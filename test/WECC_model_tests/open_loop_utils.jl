using NetworkDynamics
using OrderedCollections
using OrdinaryDiffEqRosenbrock
using DiffEqCallbacks
using CairoMakie

function _lerp(ts, vs, t)
    i = searchsortedfirst(ts, t)
    i > lastindex(ts) && return last(vs)
    i == firstindex(ts) && return first(vs)
    t0, t1 = ts[i-1], ts[i]
    v0, v1 = vs[i-1], vs[i]
    v0 + (v1 - v0) * (t - t0) / (t1 - t0)
end

"""
    voltage_source_bus(times, ur_data, ui_data; name=:vs) -> VertexModel

Vertex 2 im Open-Loop-Test: gibt eine vorgegebene Spannungstrajektorie via
linearer Interpolation aus (z.B. pvr/pvi aus PowerFactory).

Der DUT (current_source=true Bus) empfängt diese Spannung als Eingang
und berechnet daraus den Strom als Ausgang.
"""
function voltage_source_bus(times, ur_data, ui_data; name=:vs)
    ur_v = collect(Float64, ur_data)
    ui_v = collect(Float64, ui_data)
    ts_v = collect(Float64, times)
    g = (out, ins, p, t) -> begin
        out[1] = _lerp(ts_v, ur_v, t)
        out[2] = _lerp(ts_v, ui_v, t)
        nothing
    end
    VertexModel(;
        f=nothing, g, dim=0,
        outsym=[:busbar₊u_r, :busbar₊u_i],
        insym=[:busbar₊i_r, :busbar₊i_i],
        ff=NoFeedForward(), name,
    )
end

"""
    open_loop_bus_test(dut, df, voltage_rows, current_rows, observable_mapping; kwargs...)

Simuliert einen `current_source=true` Bus im Open-Loop-Modus:

- `voltage_rows`: Tupel (ur, ui) — Spaltennamen in `df` für die vorgeschriebene
  Klemmenspannung (z.B. `("pvr", "pvi")`). Wird dem DUT als Eingang vorgegeben.
- `current_rows`: Tupel (ir, ii) — Spaltennamen in `df` für den Referenzstrom
  (z.B. `("pir", "pii")`). Wird als Vergleich gegen den berechneten DUT-Ausgang geplottet.
- `observable_mapping`: `OrderedDict{Symbol,String}` — mappt interne Modell-Symbole
  auf Spalten in `obs_df` (oder `df`), die als Vergleichsreferenz angezeigt werden.
- `obs_df`: optionaler separater DataFrame für die Observable-Referenzen.

Gibt `(; sol, dat)` zurück.
"""
function open_loop_bus_test(
    dut::VertexModel, df, voltage_rows, current_rows,
    observable_mapping::AbstractDict=OrderedDict{Symbol,String}();
    solver=Rodas5P(),
    init_tol=1e-3,
    obs_df=nothing,
    external_measurements=false,
)
    cur, cui = voltage_rows
    cir, cii = current_rows
    ur_data = df[!, Symbol(cur)]
    ui_data = df[!, Symbol(cui)]
    ir_data = df[!, Symbol(cir)]
    ii_data = df[!, Symbol(cii)]
    ur0 = ur_data[1]; ui0 = ui_data[1]
    ir0 = ir_data[1]; ii0 = ii_data[1]

    # current_source=true Bus: Spannung ist Eingang, Strom ist Ausgang.
    # busbar₊u_r/u_i als Eingangs-Constraint für initialize_component!.
    # busbar₊i_r/i_i als Startwert für den Ausgang (busbar₊i_r = -pir in Netzkonvention).
    set_default!(dut, :busbar₊u_r,  ur0)
    set_default!(dut, :busbar₊u_i,  ui0)
    set_default!(dut, :busbar₊i_r, -ir0)
    set_default!(dut, :busbar₊i_i, -ii0)
    if external_measurements
        set_default!(dut, :PV₊pir_ext, -ir0)
        set_default!(dut, :PV₊pii_ext, -ii0)
    end

    residual = Ref(NaN)
    try
        initialize_component!(dut; verbose=true, tol=init_tol, residual)
    catch e
        e isa NetworkDynamics.ComponentInitError || rethrow(e)
        @warn "Initialization residuaesidual[]) above init_tol=$init_tol — proceeding anyway."
    end

    vs   = voltage_source_bus(df.time, ur_data, ui_data)
    # src=1 (DUT), dst=2 (Spannungsquelle):
    # potential von dst=2 → src=1: Spannungsquelle gibt u_r/u_i an DUT-Eingang ✓
    # flow    von src=1 → dst=2: DUT-Ausgangsstrom geht in Spannungsquelle (ignoriert) ✓
    edge = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=1, dst=2)
    nw   = Network([dut, vs], [edge])

    s0    = NWState(nw)
    tspan = (df.time[begin], df.time[end])
    prob  = ODEProblem(nw, uflat(s0), tspan, pflat(s0))

    if external_measurements
        ts_v      = collect(Float64, df.time)
        pir_ext_v = .-(ir_data)   # busbar convention: pir_model = -csv.pir
        pii_ext_v = .-(ii_data)   # busbar convention: pii_model = -csv.pii
        meas_cb = FunctionCallingCallback(
            (u, t, integrator) -> begin
                p = NWParameter(integrator)
                p.v[1, :PV₊pir_ext] = _lerp(ts_v, pir_ext_v, t)
                p.v[1, :PV₊pii_ext] = _lerp(ts_v, pii_ext_v, t)
            end;
            func_start = true,
        )
        sol = solve(prob, solver; tstops=ts_v, callback=meas_cb)
    else
        sol = solve(prob, solver; tstops=collect(df.time))
    end

    results = OrderedDict{String,Any}()
    tsim = _ol_refine(sol.t, 5)

    # Spannung: Ausgang der Spannungsquelle (trivial korrekt, zur Visualisierung)
    # Strom: Ausgang des DUT — das eigentliche Modellergebnis, verglichen mit CSV
    for (sym, col, type) in [
        (VIndex(1, :PV₊pir),     cir, 1),
        (VIndex(1, :PV₊pii),     cii, 1),
        (VIndex(2, :busbar₊u_r), cur, 2),
        (VIndex(2, :busbar₊u_i), cui, 2),
    ]
        results[col] = (;
            type, tref=df.time, ref=df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u,
            sym=string(sym.subidx),
        )
    end

    _obs_df = isnothing(obs_df) ? df : obs_df
    for (sym, col) in observable_mapping
        results[string(col)] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=VIndex(1, sym)).u,
            sym=string(sym),
        )
    end

    (; sol, dat=results)
end

"""
    comparison_figure(results; tmin=-Inf, tmax=Inf) -> Figure

Zeigt Referenz- vs. Simulationskurven für alle Einträge in `results.dat`.
Mit `tmin`/`tmax` kann der angezeigte Zeitbereich eingeschränkt werden.
"""
function comparison_figure(results; tmin=-Inf, tmax=Inf)
    res    = results.dat
    nplots = length(res)
    ncols  = max(1, floor(Int, sqrt(nplots)))
    nrows  = ceil(Int, nplots / ncols)
    fig    = Figure(size=(600 * ncols, 250 * nrows))
    idx    = 1
    for (name, r) in pairs(res)
        row = ceil(Int, idx / ncols); col = mod1(idx, ncols)
        ax  = Axis(fig[row, col]; xlabel="Time (s)", title=name)
        iref = (r.tref .>= tmin) .& (r.tref .<= tmax)
        isim = (r.tsim .>= tmin) .& (r.tsim .<= tmax)
        lines!(ax, r.tref[iref], r.ref[iref]; color=Cycled(r.type),     label="Reference",
               linewidth=2, alpha=0.7)
        lines!(ax, r.tsim[isim], r.sim[isim]; color=Cycled(r.type + 3), label=r.sym,
               linewidth=2, linestyle=:dash)
        if hasproperty(r, :extra)
            iex = (r.extra.tref .>= tmin) .& (r.extra.tref .<= tmax)
            lines!(ax, r.extra.tref[iex], r.extra.ref[iex]; color=:black,
                   label=r.extra.label, linewidth=2)
        end
        axislegend(ax; position=:lb)
        idx += 1
    end
    fig
end

function _ol_refine(ts, n)
    out = Float64[]
    for i in 1:length(ts)-1
        a, b = ts[i], ts[i+1]
        for k in 0:n-1
            push!(out, a + (b-a)*k/n)
        end
    end
    push!(out, last(ts))
    out
end
