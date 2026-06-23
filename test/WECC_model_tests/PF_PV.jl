using PowerDynamics
using OpPoDyn
using OpPoDyn.Library

using PowerDynamics.Library
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

using CSV
using DataFrames
using CairoMakie
using Test
using Printf


ref_pv = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf","variables-testcase2Bus-with-event_withoutThresholdOperationVoltage_correctSystemBase_moreTimesteps.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    drop = (i, name) -> contains(string(name), "nrows="),
    silencewarnings = true
)

ref_measures = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf","measures-testcase2Bus-with-event_staticLoad_SystemBase_moreTimesteps.csv"),
    DataFrame;
    header = 3,
    delim = ';',
    decimal = ',',
    silencewarnings = true
)
rename!(ref_measures, 1 => :time)


# bus 1 is provided from outside
PV_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.001047
    #angle_0 = 0.1
    P_0 = 0.8888
    Q_0 = -0.3333

    @named PV = OpPoDyn.Library.WECC_large_PV_pf()
    busmodel = compile_bus(MTKBus(PV); current_source=true)
    compile_bus(busmodel, pf=pfPQ(P=P_0, Q=Q_0; current_source=true)) #, assume_io_coupling=true
end

sol_pv = OpenIPSL_RePSSE_pv_pf(PV_BUS; ω_b = 2π*50);

# ── Vorfehlerzustand-Vergleich Julia vs. PowerFactory ────────────────────────
let t0 = 0.99
    i_pf = argmin(abs.(ref_pv.time .- t0))

    check = [
        (VIndex(:GEN1, :PV₊V_t),          "reec_Vt",    "V_t        "),
        (VIndex(:GEN1, :PV₊pvr),          "repc_vregr", "pvr        "),
        (VIndex(:GEN1, :PV₊pvi),          "repc_vregi", "pvi        "),
        (VIndex(:GEN1, :PV₊repca₊Q_ext),  "repc_Qext",  "Q_ext      "),
        (VIndex(:GEN1, :PV₊repca₊P_ref),  "repc_Pref",  "P_ref      "),
        (VIndex(:GEN1, :PV₊reecb₊I_qcmd), "reec_Iqcmd", "I_qcmd     "),
        (VIndex(:GEN1, :PV₊reecb₊I_pcmd), "reec_Ipcmd", "I_pcmd     "),
        (VIndex(:GEN1, :PV₊regca₊I_q),    "regc_Iq",    "regca_I_q  "),
        (VIndex(:GEN1, :PV₊regca₊I_p),    "regc_Ip",    "regca_I_p  "),
        (VIndex(:GEN1, :PV₊regca₊I_qr),   "regc_Iqr",   "regca_I_qr "),
    ]

    println("\n=== Vorfehlerzustand t=$(t0)s: Julia vs. PowerFactory ===")
    @printf("%-16s  %12s  %12s  %12s\n", "Variable", "Julia", "PF", "Diff")
    for (idx, col, name) in check
        jv  = sol_pv(t0, idxs=idx)
        col_sym = Symbol(col)
        if col_sym ∉ propertynames(ref_pv)
            @printf("%-16s  %12.6f  %12s  %12s\n", name, jv, "(no col)", "—")
            continue
        end
        pv = ref_pv[i_pf, col_sym]
        @printf("%-16s  %12.6f  %12.6f  %+.2e\n", name, jv, pv, jv - pv)
    end
end

# ── Post-Fault-Vergleich Julia vs. PowerFactory ───────────────────────────
# t=1.001s: Strom und V_t direkt nach dem Fault (algebraische Reinitialisierung)
let t0 = 1.001
    i_pf = argmin(abs.(ref_pv.time .- t0))
    check = [
        (VIndex(:GEN1, :PV₊V_t),          "reec_Vt",      "V_t     "),
        (VIndex(:GEN1, :PV₊Q_measure),    "repc_Qbranch", "Q_branch"),
        (VIndex(:GEN1, :PV₊P_measure),    "repc_Pbranch", "P_branch"),
        (VIndex(:GEN1, :PV₊reecb₊I_qcmd), "reec_Iqcmd",   "I_qcmd  "),
        (VIndex(:GEN1, :PV₊reecb₊I_pcmd), "reec_Ipcmd",   "I_pcmd  "),
        (VIndex(:GEN1, :PV₊regca₊I_q),    "regc_Iq",      "regca_Iq"),
        (VIndex(:GEN1, :PV₊regca₊I_p),    "regc_Ip",      "regca_Ip"),
        (VIndex(:GEN1, :PV₊regca₊I_qr),   "regc_Iqr",     "regca_Iqr"),
    ]
    println("\n=== Post-Fault t=$(t0)s (Stromzustände + V_t) ===")
    @printf("%-12s  %12s  %12s  %12s\n", "Variable", "Julia", "PF", "Diff")
    for (idx, col, name) in check
        jv = sol_pv(t0, idxs=idx)
        col_sym = Symbol(col)
        col_sym ∉ propertynames(ref_pv) && (@printf("%-12s  %12.6f  %12s\n", name, jv, "(no col)"); continue)
        pv = ref_pv[i_pf, col_sym]
        @printf("%-12s  %12.6f  %12.6f  %+.2e\n", name, jv, pv, jv - pv)
    end
end

# t=1.02s und t=1.3s: zeitliche Entwicklung der Abweichung
let check = [
        (VIndex(:GEN1, :PV₊V_t),          "reec_Vt",      "V_t     "),
        (VIndex(:GEN1, :PV₊Q_measure),    "repc_Qbranch", "Q_branch"),
        (VIndex(:GEN1, :PV₊P_measure),    "repc_Pbranch", "P_branch"),
        (VIndex(:GEN1, :PV₊reecb₊I_pcmd), "reec_Ipcmd",   "I_pcmd  "),
        (VIndex(:GEN1, :PV₊regca₊I_p),    "regc_Ip",      "regca_Ip"),
    ]
    for t0 in [1.02, 1.3]
        i_pf = argmin(abs.(ref_pv.time .- t0))
        println("\n=== Post-Fault t=$(t0)s ===")
        @printf("%-10s  %12s  %12s  %12s\n", "Variable", "Julia", "PF", "Diff")
        for (idx, col, name) in check
            jv = sol_pv(t0, idxs=idx)
            col_sym = Symbol(col)
            col_sym ∉ propertynames(ref_pv) && continue
            pv = ref_pv[i_pf, col_sym]
            @printf("%-10s  %12.6f  %12.6f  %+.2e\n", name, jv, pv, jv - pv)
        end
    end
end

## perform tests for all variables of interest
# Plant controls (repc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊P_ref), "repc_Pref") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊Q_ext), "repc_Qext") < 1e-3

# Electrical control (reec_b)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊Q_measure), "repc_Qbranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊P_measure), "repc_Pbranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊V_t), "reec_Vt") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pcmd), "reec_Ipcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qcmd), "reec_Iqcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pmax), "reec_Ipmax0") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmax), "reec_Iqmax0") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmin), "reec_Iqmin0") < 1e-3

# Renewable generator (regc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊regca₊I_lvpl), "regc_Ilvpl") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊I_measure), "repc_Ibranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvi), "repc_vregi") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvr), "repc_vregr") < 1e-3

#Input Plot
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig2 = let
        fig = Figure(resolution=(1400, 1200))
        ts   = range(0.0, 1.5; length=2000)
        xlims = (0.0, 1.5)

        # Plot 1: pvr
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="pvr (V real)", limits=(xlims..., nothing, nothing))
        lines!(ax1, ref_measures.time, ref_measures.ur; label="PowerFactory", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="Julia", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: pvi
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="pvi (V imag)", limits=(xlims..., nothing, nothing))
        lines!(ax2, ref_measures.time, ref_measures.ui; label="PowerFactory", color=:green, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="Julia", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: V_t
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="V_t (terminal voltage magnitude)", limits=(xlims..., nothing, nothing))
        lines!(ax3, ref_measures.time, ref_measures.u; label="PowerFactory", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="Julia", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: I_branch
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="I_branch (current magnitude)", limits=(xlims..., nothing, nothing))
        lines!(ax4, ref_measures.time, ref_measures.i; label="PowerFactory", color=:red, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊I_measure)).u; label="Julia", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        # Plot 5: P_branch
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="P_branch (active power)", limits=(xlims..., nothing, nothing))
        lines!(ax5, ref_measures.time, ref_measures.p; label="PowerFactory", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊P_measure)).u; label="Julia", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        # Plot 6: Q_branch
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="[pu]", title="Q_branch (reactive power)", limits=(xlims..., nothing, nothing))
        lines!(ax6, ref_measures.time, ref_measures.q; label="PowerFactory", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u; label="Julia", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        # Plot 7: pir (terminal current real part)
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="[pu]", title="pir (I real)", limits=(xlims..., nothing, nothing))
        lines!(ax7, ref_measures.time, ref_measures.ir; label="PowerFactory", color=:red, linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pir)).u; label="Julia", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax7)

        # Plot 8: pii (terminal current imaginary part)
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="[pu]", title="pii (I imag)", limits=(xlims..., nothing, nothing))
        lines!(ax8, ref_measures.time, ref_measures.ii; label="PowerFactory", color=:teal, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pii)).u; label="Julia", color=:teal, linestyle=:dash, linewidth=2)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_comparison_network.png"), fig2)
end

#plot output submodules
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig3 = let
        fig = Figure(resolution=(1400, 1000))
        ts    = range(0.0, 1.5; length=2000)
        xlims = (0.0, 1.5)

        # Plot 1: P_ref
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="P_ref", limits=(xlims..., nothing, nothing))
        lines!(ax1, ref_pv.time, ref_pv[!, Symbol("repc_Pref")]; label="PowerFactory", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊P_ref)).u; label="Julia", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: Q_ext
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="Q_ext", limits=(xlims..., nothing, nothing))
        lines!(ax2, ref_pv.time, ref_pv[!, Symbol("repc_Qext")]; label="PowerFactory", color=:red, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_ext)).u; label="Julia", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: I_pcmd
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="I_pcmd", limits=(xlims..., nothing, nothing))
        lines!(ax3, ref_pv.time, ref_pv[!, Symbol("reec_Ipcmd")]; label="PowerFactory", color=:green, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_pcmd)).u; label="Julia", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: I_qcmd
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="I_qcmd", limits=(xlims..., nothing, nothing))
        lines!(ax4, ref_pv.time, ref_pv[!, Symbol("reec_Iqcmd")]; label="PowerFactory", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_qcmd)).u; label="Julia", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        # Plot 5: I_pout
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="[pu]", title="I_pout", limits=(xlims..., nothing, nothing))
        lines!(ax5, ref_pv.time, ref_pv[!, Symbol("regc_Ip")]; label="PowerFactory", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_p)).u; label="Julia", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        # Plot 6: I_qout
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="[pu]", title="I_qout", limits=(xlims..., nothing, nothing))
        lines!(ax6, ref_pv.time, ref_pv[!, Symbol("regc_Iq")]; label="PowerFactory", color=:teal, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_q)).u; label="Julia", color=:teal, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_comparison_controls.png"), fig3)
end


# ══════════════════════════════════════════════════════════════════════════════
# Open-Loop-Test: Terminal-Größen aus measures-CSV vorgeben, interne Zustände prüfen
# pvr/pvi/pir/pii aus PowerFactory → WECC-Modell → interne Größen vergleichen
# ══════════════════════════════════════════════════════════════════════════════
include(joinpath(@__DIR__, "open_loop_utils.jl"))
include(joinpath(@__DIR__, "prescribed_pv_test_2bus.jl"))

# Spaltennamen in das von prescribed_pv_test_2bus erwartete Format umbenennen
ref_inputs_2bus = let
    df = copy(ref_measures)
    rename!(df,
        :ur    => :pvr,
        :ui    => :pvi,
        :ir    => :pir,
        :ii    => :pii,
        :i     => :I_measure,
        :p     => :P_measure,
        :q     => :Q_measure,
        :u     => :V_t,
        :fmeas => :f,
    )
    df
end

res_presc = prescribed_pv_test_2bus(
    ref_inputs_2bus;
    obs_df = ref_pv,
)

regc_keys        = ["regc_Iqr", "regc_o2", "regc_o3", "regc_dIq", "regc_Qgen", "regc_Vtfiltlim"]
repc_q_diag_keys = ["repc_Qfltr", "repc_dQ", "repc_dQin", "repc_dQdbd", "repc_Qe", "repc_Qlim"]
exclude_keys     = union(regc_keys, repc_q_diag_keys)

fig_presc = comparison_figure(
    (; sol=res_presc.sol, dat=OrderedDict(k => v for (k, v) in res_presc.dat if k ∉ exclude_keys));
    tmin=0.0, tmax=1.5,
)

fig_regc_diag = comparison_figure(
    (; sol=res_presc.sol, dat=OrderedDict(k => res_presc.dat[k] for k in regc_keys));
    tmin=0.0, tmax=1.5,
)

repc_q_keys = ["repc_Qbranch", "repc_Qfltr", "repc_dQ", "repc_dQin",
               "repc_dQdbd", "repc_Qe", "repc_Qlim", "repc_Qext"]
fig_repc_q_diag = comparison_figure(
    (; sol=res_presc.sol, dat=OrderedDict(k => res_presc.dat[k] for k in repc_q_keys));
    tmin=0.0, tmax=2.0,
)
save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_2bus_repc_q_diag.png"), fig_repc_q_diag)

# Focused Q-path plot: Q_branch → ΔQ_in → ΔQ_dbd → Q_e → Q_lim → Q_ext (3×2)
fig_repc_q_focus = let
    keys6  = ["repc_Qbranch", "repc_dQin", "repc_dQdbd", "repc_Qe", "repc_Qlim", "repc_Qext"]
    titles6 = ["Q_branch (measured)", "ΔQ_in (after deadband input)", "ΔQ_dbd (deadband output)",
               "Q_e (error signal)", "Q_lim (PI output)", "Q_ext (REPC output)"]
    fig = Figure(size=(1100, 900))
    for (i, (k, title)) in enumerate(zip(keys6, titles6))
        r   = res_presc.dat[k]
        row = div(i - 1, 2) + 1
        col = mod(i - 1, 2) + 1
        ax  = Axis(fig[row, col]; xlabel="Time (s)", ylabel="[pu]", title=title)
        iref = (r.tref .>= 0.0) .& (r.tref .<= 2.0)
        isim = (r.tsim .>= 0.0) .& (r.tsim .<= 2.0)
        lines!(ax, r.tref[iref], r.ref[iref]; color=:blue,   label="PowerFactory", linewidth=2, linestyle=:solid)
        lines!(ax, r.tsim[isim], r.sim[isim]; color=:orange, label="Julia",        linewidth=2)
        axislegend(ax; position=:lt)
    end
    fig
end

# Fehlerplots: Differenz Julia − PowerFactory für Q_lim und ΔQ_in
fig_qlim_error = let
    fig = Figure(size=(900, 600))

    for (row, (key, ylabel, title)) in enumerate([
        ("repc_Qlim", "Δ Q_lim [pu]",  "Q_lim: Julia − PowerFactory"),
        ("repc_dQin", "Δ ΔQ_in [pu]",  "ΔQ_in: Julia − PowerFactory"),
    ])
        r    = res_presc.dat[key]
        isim = (r.tsim .>= 0.0) .& (r.tsim .<= 2.5)
        t    = r.tsim[isim]
        diff = r.sim[isim] .- _lerp.(Ref(r.tref), Ref(r.ref), t)

        ax = Axis(fig[row, 1]; xlabel="Time (s)", ylabel=ylabel, title=title)
        hlines!(ax, 0.0; color=:black, linewidth=1, linestyle=:dot)
        lines!(ax, t, diff; color=:red, linewidth=2)
    end

    fig
end

# Q-Pfad Upstream: Q_branch, Q_fltr, ΔQ, ΔQ_in — Julia vs. PowerFactory
fig_qpath_upstream = let
    entries = [
        ("repc_Qbranch", "Q_branch [pu]",  "Q_branch (gemessen)"),
        ("repc_Qfltr",   "Q_fltr [pu]",    "Q_fltr (nach SimpleLag)"),
        ("repc_dQ",      "ΔQ [pu]",        "ΔQ = Q_ref − Q_fltr"),
        ("repc_dQin",    "ΔQ_in [pu]",     "ΔQ_in (nach RefFlag)"),
    ]
    fig = Figure(size=(1100, 900))
    for (i, (key, ylabel, title)) in enumerate(entries)
        r    = res_presc.dat[key]
        row  = div(i - 1, 2) + 1
        col  = mod(i - 1, 2) + 1
        iref = (r.tref .>= 0.0) .& (r.tref .<= 2.5)
        isim = (r.tsim .>= 0.0) .& (r.tsim .<= 2.5)
        ax   = Axis(fig[row, col]; xlabel="Time (s)", ylabel=ylabel, title=title)
        lines!(ax, r.tref[iref], r.ref[iref]; color=:blue,   label="PowerFactory", linewidth=2)
        lines!(ax, r.tsim[isim], r.sim[isim]; color=:orange, label="Julia",        linewidth=2)
        axislegend(ax; position=:lt)
    end
    fig
end

# Fehlerplots: Q_branch und Q_fltr — Julia − PowerFactory
fig_qbranch_fltr_error = let
    fig = Figure(size=(900, 600))
    for (row, (key, ylabel, title)) in enumerate([
        ("repc_Qbranch", "Δ Q_branch [pu]", "Q_branch: Julia − PowerFactory"),
        ("repc_Qfltr",   "Δ Q_fltr [pu]",   "Q_fltr: Julia − PowerFactory"),
    ])
        r    = res_presc.dat[key]
        isim = (r.tsim .>= 0.0) .& (r.tsim .<= 2.5)
        t    = r.tsim[isim]
        diff = r.sim[isim] .- _lerp.(Ref(r.tref), Ref(r.ref), t)
        ax   = Axis(fig[row, 1]; xlabel="Time (s)", ylabel=ylabel, title=title)
        hlines!(ax, 0.0; color=:black, linewidth=1, linestyle=:dot)
        lines!(ax, t, diff; color=:red, linewidth=2)
    end
    fig
end
