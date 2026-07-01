# Standalone open-loop test for WECC_large_PV_prescribed.
#
# All 4 terminal quantities (pvr, pvi, pir, pii) are prescribed from CSV via
# registered signal functions (_wecc_pvr_sig etc.).  No callbacks needed —
# the functions simply interpolate the stored data at whatever t the solver
# asks for.  Call set_prescribed_terminal! once before building the problem.
#
# Sign convention: pir/pii taken directly from CSV (same sign as terminal.i_r in the MTK model).
# No negation needed here — the @register_symbolic signals feed the MTK model directly,
# bypassing the NetworkDynamics busbar layer where the sign flip (busbar.i_r = -csv.pir) lives.

using SciMLBase: SciMLBase

function prescribed_pv_test(
    df;
    solver         = Rodas5P(autodiff=false),
    obs_df         = nothing,
    calc_measures  = false,   # if true: I/P/Q_measure computed from pvr/pvi/pir/pii instead of CSV
    P0             = 0.8888,  # initial active power reference (overrides model default)
    Q0             = -0.3333, # initial reactive power reference (overrides model default)
)
    pvr_d   = collect(Float64, df[!, :pvr])
    pvi_d   = collect(Float64, df[!, :pvi])
    pir_d   = collect(Float64, df[!, :pir])   # same sign as terminal.i_r in the MTK model
    pii_d   = collect(Float64, df[!, :pii])
    imeas_d = collect(Float64, df[!, :I_measure])
    pmeas_d = collect(Float64, df[!, :P_measure])
    qmeas_d = collect(Float64, df[!, :Q_measure])
    ts      = collect(Float64, df.time)

    # Load CSV data into the Library-level signal store
    OpPoDyn.Library.set_prescribed_terminal!(ts, pvr_d, pvi_d, pir_d, pii_d, imeas_d, pmeas_d, qmeas_d)

    @named PV = OpPoDyn.Library.WECC_large_PV_prescribed(; calc_measures)
    sys   = structural_simplify(PV)
    tspan = (ts[1], ts[end])

    # Build explicit u0 from the model's guess values.  MTK's initialization is
    # underdetermined for @register_symbolic inputs (opaque to structural analysis),
    # so we skip it entirely with NoInit() and rely on the guess values instead.
    g  = ModelingToolkit.get_guesses(sys)
    u0 = Float64[get(g, s, 0.0) for s in unknowns(sys)]

    # Compute initial terminal voltage magnitude from CSV
    V0 = sqrt(pvr_d[1]^2 + pvi_d[1]^2)

    # Override initial states so the model starts in steady state at (P0, Q0)
    # instead of the hardcoded 2-bus guess values (0.800721, -0.30027)
    OLD_Q0 = -0.30027
    for (i, s) in enumerate(unknowns(sys))
        sstr = string(s)
        # REPC: Q_fltr, PI integral, leadLag → Q0
        if contains(sstr, "repca") &&
           (contains(sstr, "simpleLag2") || contains(sstr, "PI_lim_Q") || contains(sstr, "leadLag"))
            if get(g, s, 0.0) ≈ OLD_Q0
                u0[i] = Q0
            end
        # REEC: P_refout lag (P_limLag) → P0
        elseif contains(sstr, "reecb") && contains(sstr, "P_limLag")
            u0[i] = P0
        # REEC: I_qin lag (simpleLag_freeze) → Q0/V0
        elseif contains(sstr, "reecb") && contains(sstr, "simpleLag_freeze")
            u0[i] = Q0/V0
        # REGCA: I_pr lag with rate limit (SimpleLag_2uplims) → P0/V0
        elseif contains(sstr, "regca") && contains(sstr, "SimpleLag_2uplims")
            u0[i] = P0/V0
        # REGCA: I_qr lag (SimpleLagLim) → Q0/V0
        elseif contains(sstr, "regca") && contains(sstr, "SimpleLagLim")
            u0[i] = Q0/V0
        end
    end

    # _callback_sat_min/_max: internal limiter parameters not set by default
    sat_defaults = [p => 0.0 for p in parameters(sys)
                    if contains(string(p), "_callback_sat")]
    # Override plant-level initial condition parameters (P0/Q0)
    init_params = [p => (contains(string(p), "Qref_set") || contains(string(p), "Qinit_set") ? Q0 : P0)
                   for p in parameters(sys)
                   if any(s -> contains(string(p), s),
                          ["P_plantref_set", "Pinit_set", "Qref_set", "Qinit_set"])]
    prob = ODEProblem(sys, u0, tspan, vcat(sat_defaults, init_params))
    prob = PowerDynamics.Library.attach_limint_callbacks_ode(sys, prob)

    sol  = solve(prob, solver; tstops=ts, initializealg=SciMLBase.NoInit())
    @info "prescribed_pv_test: retcode=$(sol.retcode), nsteps=$(length(sol.t))"
    tsim = _ol_refine(sol.t, 5)

    results = OrderedDict{String,Any}()
    _obs_df = isnothing(obs_df) ? df : obs_df

    # ── 1. Terminal inputs ────────────────────────────────────────────────────
    V_t_sim = sol(tsim; idxs=sys.V_t).u
    δ_v_sim = sol(tsim; idxs=sys.δ_v).u
    results["pvr"] = (;
        type=2, tref=ts, ref=pvr_d,
        tsim, sim=V_t_sim .* cos.(δ_v_sim), sym="V_t·cos(δ_v)")
    results["pvi"] = (;
        type=2, tref=ts, ref=pvi_d,
        tsim, sim=V_t_sim .* sin.(δ_v_sim), sym="V_t·sin(δ_v)")
    results["pii"] = (;
        type=1, tref=ts, ref=pii_d,
        tsim, sim=_lerp.(Ref(ts), Ref(pii_d), tsim), sym="pii (signal)")
    results["pir"] = (;
        type=1, tref=ts, ref=pir_d,
        tsim, sim=_lerp.(Ref(ts), Ref(pir_d), tsim), sym="pir (signal)")

    # ── 2. Measurements (directly from ref_inputs — same source as repc_Ibranch/Qbranch/Pbranch) ──
    I_measure_d = collect(Float64, df[!, :I_measure])
    P_measure_d = collect(Float64, df[!, :P_measure])
    Q_measure_d = collect(Float64, df[!, :Q_measure])
    results["repc_Ibranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Ibranch],
        tsim, sim=sol(tsim; idxs=sys.I_measure).u, sym="I_measure")
    results["repc_Pbranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Pbranch],
        tsim, sim=sol(tsim; idxs=sys.P_measure).u, sym="P_measure")
    results["repc_Qbranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Qbranch],
        tsim, sim=sol(tsim; idxs=sys.Q_measure).u, sym="Q_measure")

    # ── 3. REPC ───────────────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.repca.Q_ext, "repc_Qext"),
        (sys.repca.P_ref, "repc_Pref"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 4. REEC ───────────────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.reecb.I_qcmd, "reec_Iqcmd"),
        (sys.reecb.I_pcmd, "reec_Ipcmd"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 5. REGC ───────────────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.regca.I_q, "regc_Iq"),
        (sys.regca.I_p, "regc_Ip"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 6. Model output current ───────────────────────────────────────────────
    results["pii_out"] = (;
        type=1, tref=ts, ref=pii_d,
        tsim, sim=sol(tsim; idxs=sys.pii_out).u, sym="pii_out")
    results["pir_out"] = (;
        type=1, tref=ts, ref=pir_d,
        tsim, sim=sol(tsim; idxs=sys.pir_out).u, sym="pir_out")

    # ── 7. REGC internals ────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.regca.I_qr,       "regc_Iqr"),
        (sys.regca.o2,         "regc_o2"),
        (sys.regca.o3,         "regc_o3"),
        (sys.regca.ΔI_q,       "regc_dIq"),
        (sys.regca.Q_gen,      "regc_Qgen"),
        (sys.regca.V_tfiltlim, "regc_Vtfiltlim"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 8. REPC Q_ext Signalpfad ─────────────────────────────────────────────
    for (sym, col) in [
        (sys.repca.Q_fltr,  "repc_Qfltr"),
        (sys.repca.ΔQ,      "repc_dQ"),
        (sys.repca.ΔQ_in,   "repc_dQin"),
        (sys.repca.ΔQ_dbd,  "repc_dQdbd"),
        (sys.repca.Q_e,     "repc_Qe"),
        (sys.repca.Q_lim,   "repc_Qlim"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    (; sol, dat=results, sys)
end
