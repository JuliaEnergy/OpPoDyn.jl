# Standalone open-loop test for the 2-bus WECC PV case.
#
# Uses WECC_large_PV_prescribed with 2-bus-specific parameter and initial-condition
# overrides.  After structural_simplify, sub-system property access (sys.reecb.X) can
# fail for states that were reorganised by MTK.  We therefore use string matching on
# unknowns(sys) — more robust across MTK versions.
#
# State-name conventions after structural_simplify (MTK ₊ separator):
#   PowerDynamics.SimpleLag        → "…₊internal(t)"
#   PowerDynamics.LeadLag          → "…₊internal(t)"
#   PowerDynamics.SimpleLagLim     → "…₊x(t)"
#   OpPoDyn.P_I_Lim_freeze         → "…₊x(t)"
#   OpPoDyn.SimpleLag_2MaxLims     → "…₊x(t)"
#   OpPoDyn.SimpleLag_freeze       → "…₊x(t)"
#   OpPoDyn.SimpleLag_2Lims_freeze → "…₊x(t)"
#
# reec_b_pf specifics for QFlag=false:
#   - I_qin is NOT a state; it equals simpleLag_freeze.out (computed)
#   - P_refout is commented out; P lag is P_limLag (SimpleLag_2Lims_freeze)
#   - s_Qint / PI_freeze do NOT exist for QFlag=false (Q→V block absent)

using SciMLBase: SciMLBase

# ── string-matching override helper ─────────────────────────────────────────────
# Sets u0[i] = val for every state whose string representation contains ALL
# patterns in `include` and NONE of the patterns in `exclude`.
function _override_match!(u0, sys, include::Vector{String}, val;
                          exclude::Vector{String}=String[])
    sts  = unknowns(sys)
    strs = string.(sts)
    idxs = findall(i ->
        all(p -> occursin(p, strs[i]), include) &&
        all(p -> !occursin(p, strs[i]), exclude),
        eachindex(sts))
    if isempty(idxs)
        @warn "prescribed_pv_test_2bus: no state matching $include (excl $exclude) — IC not set"
    else
        u0[idxs] .= val
        @info "prescribed_pv_test_2bus: set $(strs[idxs]) = $val"
    end
end


function prescribed_pv_test_2bus(
    df;
    solver         = Rodas5P(autodiff=false),
    obs_df         = nothing,
    calc_measures  = false,
    Vref_set       = 0.903563,
    Qref_set       = -0.3333,
    P_plantref_set = 0.8888,
    Qinit_set      = -0.3333,
    Pinit_set      = 0.8888,
)
    pvr_d   = collect(Float64, df[!, :pvr])
    pvi_d   = collect(Float64, df[!, :pvi])
    pir_d   = collect(Float64, df[!, :pir])
    pii_d   = collect(Float64, df[!, :pii])
    imeas_d = collect(Float64, df[!, :I_measure])
    pmeas_d = collect(Float64, df[!, :P_measure])
    qmeas_d = collect(Float64, df[!, :Q_measure])
    ts      = collect(Float64, df.time)

    OpPoDyn.Library.set_prescribed_terminal!(ts, pvr_d, pvi_d, pir_d, pii_d, imeas_d, pmeas_d, qmeas_d)

    @named PV = OpPoDyn.Library.WECC_large_PV_prescribed(; calc_measures)
    sys = structural_simplify(PV)

    # ── 2-bus parameter overrides ─────────────────────────────────────────────
    sat_defaults = [p => 0.0 for p in parameters(sys)
                    if contains(string(p), "_callback_sat")]
    param_2bus = vcat(sat_defaults, [
        sys.Vref_set       => Vref_set,
        sys.Qref_set       => Qref_set,
        sys.P_plantref_set => P_plantref_set,
        sys.Qinit_set      => Qinit_set,
        sys.Pinit_set      => Pinit_set,
    ])

    # ── Steady-state values for 2-bus operating point ─────────────────────────
    Q_0   = Float64(Qref_set)        # -0.3333
    P_0   = Float64(P_plantref_set)  # 0.8888
    V_t_0 = Float64(Vref_set)        # 0.903563
    I_q0  = Q_0 / V_t_0             # ≈ -0.3690  (I_qin and I_q filter steady-state)
    I_p0  = P_0 / V_t_0             # ≈ 0.9836   (I_p filter steady-state)

    # ── Initial conditions: guesses + targeted 2-bus overrides ───────────────
    g  = ModelingToolkit.get_guesses(sys)
    u0 = Float64[get(g, s, 0.0) for s in unknowns(sys)]

    # ── REPC (repca) ──────────────────────────────────────────────────────────
    # PI_lim_Q.x   OpPoDyn P_I_Lim_freeze, state "x"
    # FROZEN — cannot self-correct when Q-error = 0
    _override_match!(u0, sys, ["PI_lim_Q", "x"],        Q_0)
    # leadLag.internal  PowerDynamics LeadLag, state "internal"
    # Q_ext lag (T_fv=0.1s): converges to PI output → must be set explicitly
    _override_match!(u0, sys, ["leadLag", "internal"],   Q_0)
    # simpleLag2.internal  PowerDynamics SimpleLag, state "internal" — Q filter
    _override_match!(u0, sys, ["simpleLag2", "internal"], Q_0)
    # simpleLag3.internal  PowerDynamics SimpleLag, state "internal" — V filter
    _override_match!(u0, sys, ["simpleLag3", "internal"], V_t_0)

    # ── REEC (reecb) — reec_b_pf with QFlag=false ────────────────────────────────
    # simpleLag.internal (V filter, T_rv) — PowerDynamics SimpleLag, state "internal"
    # Exclude simpleLag1 (P filter); simpleLag_freeze has state "x" so no collision
    _override_match!(u0, sys, ["reecb", "simpleLag", "internal"], V_t_0;
                     exclude=["simpleLag1"])
    # simpleLag1.internal (P filter, T_p) — PowerDynamics SimpleLag, state "internal"
    _override_match!(u0, sys, ["simpleLag1", "internal"], P_0)
    # simpleLag_freeze.x (I_q lag, T_iq) — OpPoDyn SimpleLag_freeze, state "x"
    # I_qin = simpleLag_freeze.out = (Q_ext/V_t) filtered; at SS = I_q0
    _override_match!(u0, sys, ["simpleLag_freeze", "x"], I_q0)
    # P_limLag.x (P order lag, T_pord) — OpPoDyn SimpleLag_2Lims_freeze, state "x"
    # P_refout is commented out in reec_b_pf; P_lim = P_limLag.out = P_0 at SS
    _override_match!(u0, sys, ["P_limLag", "x"], P_0)
    # NOTE: for QFlag=false, the entire Q→V PI block (s_Qint, PI_freeze, etc.)
    # does not exist in reec_b_pf — reactive current = I_t_filtered directly.

    # ── REGC (regca) ──────────────────────────────────────────────────────────
    # simpleLag.internal (V filter, T_fltr=0.02s)
    _override_match!(u0, sys, ["regca", "simpleLag", "internal"], V_t_0)
    # SimpleLagLim.x  PowerDynamics LimitedIntegratorBase, state "x" — I_q lag
    _override_match!(u0, sys, ["SimpleLagLim", "x"],         I_q0)
    # SimpleLag_2uplims.x  OpPoDyn SimpleLag_2MaxLims, state "x" — I_p lag
    _override_match!(u0, sys, ["SimpleLag_2uplims", "x"],    I_p0)

    # ── Simulation ────────────────────────────────────────────────────────────
    tspan = (ts[1], ts[end])
    prob  = ODEProblem(sys, u0, tspan, param_2bus)
    sol   = solve(prob, solver; tstops=ts, initializealg=SciMLBase.NoInit())
    @info "prescribed_pv_test_2bus: retcode=$(sol.retcode), nsteps=$(length(sol.t))"
    tsim  = _ol_refine(sol.t, 5)

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

    # ── 2. Measurements ──────────────────────────────────────────────────────
    results["repc_Ibranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Ibranch],
        tsim, sim=sol(tsim; idxs=sys.I_measure).u, sym="I_measure")
    results["repc_Pbranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Pbranch],
        tsim, sim=sol(tsim; idxs=sys.P_measure).u, sym="P_measure")
    results["repc_Qbranch"] = (;
        type=3, tref=_obs_df.time, ref=_obs_df[!, :repc_Qbranch],
        tsim, sim=sol(tsim; idxs=sys.Q_measure).u, sym="Q_measure")

    # ── 3. REPC ──────────────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.repca.Q_ext, "repc_Qext"),
        (sys.repca.P_ref, "repc_Pref"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 4. REEC ──────────────────────────────────────────────────────────────
    for (sym, col) in [
        (sys.reecb.I_qcmd, "reec_Iqcmd"),
        (sys.reecb.I_pcmd, "reec_Ipcmd"),
    ]
        results[col] = (;
            type=3, tref=_obs_df.time, ref=_obs_df[!, Symbol(col)],
            tsim, sim=sol(tsim; idxs=sym).u, sym=col)
    end

    # ── 5. REGC ──────────────────────────────────────────────────────────────
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

    # ── 8. REPC Q-path diagnostics ───────────────────────────────────────────
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
