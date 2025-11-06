function OpenIPSL_RePSSE(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
    # copy constructor and set vidxs
    bus1 = VertexModel(_bus1, vidx=1, name=:GEN1)

    S_b = 100e6
    #ω_b = 2π*50 #50 bei BESS und WT4B; 60 bei PV

    bus3 = let
        # OpenIPSL infinite bus parameters from SMIB base class
        H = 0                    # H=0 makes it behave like infinite bus
        M_b = 100e6
        X_d = 0.2               # Internal impedance
        D = 0
        # V_b = 400e3

        # pf results, just used for pf modek
        # P_0 = 10.017110e6       # From OpenIPSL SMIB.mo
        # Q_0 = 8.006544e6        # From OpenIPSL SMIB.mo
        v_0 = 1.0
        angle_0 = 0 #-0.0000157           # From OpenIPSL SMIB.mo

        @named gencls_inf = PSSE_GENCLS(; S_b, ω_b, H, M_b, X_d, D)
        busmodel = MTKBus(gencls_inf; name=:GEN2)
        compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0), vidx=3)
    end

    bus2 = let
        @named pwFault = ConstantYLoad(B=0, G=0, allow_zero_conductance=true)
        busmodel = MTKBus(pwFault; name=:FAULT)
        #v_0 = 1.0001
        #angle_0 = deg2rad(0.014)
        faultbus = compile_bus(busmodel, vidx=2) #, pf=pfSlack(V=v_0, δ=angle_0)

        enable = ComponentAffect([], [:pwFault₊B, :pwFault₊G]) do u, p, ctx
            p[:pwFault₊B] = -1 #B=-X/(X^2+R^2)
            p[:pwFault₊G] = 1 #G=R/(X^2+R^2)
        end
        disable = ComponentAffect([], [:pwFault₊B, :pwFault₊G]) do u, p, ctx
            p[:pwFault₊B] = 0
            p[:pwFault₊G] = 0
        end
        enable_cb = PresetTimeComponentCallback(2, enable)
        disable_cb = PresetTimeComponentCallback(2.15, disable)
        set_callback!(faultbus, (enable_cb, disable_cb))
        faultbus
    end


    # line template
    pwLine = MTKLine(PiLine(; name=:PwLine))
    line2 = compile_line(pwLine; name=:pwLine2,
        src=:GEN1, dst=:FAULT,
        PwLine₊X=0.0025, PwLine₊R=0.0025)

    @named branchA = PiLine(; name=:pwLine,
        X=0.025, R=0.025, B_src=0.025, B_dst=0.025)
    @named branchB = PiLine(; name=:pwLine1,
        X=0.025, R=0.025, B_src=0.025, B_dst=0.025)
    linemodel = MTKLine(branchA, branchB)
    parallelline = compile_line(linemodel, src=:FAULT, dst=:GEN2)


    buses = [bus1, bus2, bus3]
    lines = [parallelline, line2]
    nw = Network(buses, lines; warn_order=false)

    verbose = true
    pfnw=nothing
    pfs0=nothing
    pfs=nothing
    pfnw = isnothing(pfnw) ? powerflow_model(nw) : pfnw
    pfs0 = isnothing(pfs0) ? NWState(pfnw) : pfnw
    pfs = solve_powerflow(nw; pfnw, pfs0, verbose)
    println(show_powerflow(pfs))
    interface_vals = interface_values(pfs)
    println(interface_vals)
    # pfnw = powerflow_model(nw)
    # pfs = solve_powerflow(pfnw)

    if just_init
        s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol=Inf, nwtol=Inf)
        return s0
    end

    for sym in sym(bus1)
        has_guess(bus1, sym) || continue
        (sym==:PV₊plant_control₊p_0) && continue
        (sym==:PV₊plant_control₊Voltage_dip) && continue
        (sym==:PV₊plant_control₊V_droop) && continue
        (sym==:PV₊plant_control₊V_in) && continue
        (sym==:PV₊plant_control₊V_fltr) && continue
        (sym==:PV₊plant_control₊ΔV) && continue
        (sym==:PV₊plant_control₊Q_fltr) && continue
        (sym==:PV₊plant_control₊ΔQ) && continue
        (sym==:PV₊plant_control₊ΔQ_in) && continue
        (sym==:PV₊plant_control₊ΔQ_dbd) && continue
        (sym==:PV₊plant_control₊Q_e) && continue
        (sym==:PV₊plant_control₊Q_x) && continue
        (sym==:PV₊plant_control₊Q_res) && continue
        (sym==:PV₊plant_control₊Q_I) && continue
        (sym==:PV₊plant_control₊Q_lim) && continue
        (sym==:PV₊plant_control₊Q_ext) && continue
        (sym==:PV₊plant_control₊Δf_deadband) && continue
        (sym==:PV₊plant_control₊Δf_corr) && continue
        (sym==:PV₊plant_control₊P_branchp) && continue
        (sym==:PV₊plant_control₊f_e) && continue
        (sym==:PV₊plant_control₊P_e) && continue
        (sym==:PV₊plant_control₊P_lim) && continue
        (sym==:PV₊plant_control₊P_refa) && continue
        (sym==:PV₊plant_control₊P_ref) && continue
        (sym==:PV₊electrical_control₊Voltage_dip) && continue
        (sym==:PV₊electrical_control₊V_tfilt) && continue
        (sym==:PV₊electrical_control₊V_tfiltlim) && continue
        (sym==:PV₊electrical_control₊ΔV_t) && continue
        (sym==:PV₊electrical_control₊ΔV_tdbd) && continue
        (sym==:PV₊electrical_control₊I_qinj) && continue
        (sym==:PV₊electrical_control₊P_PF) && continue
        (sym==:PV₊electrical_control₊Q_con) && continue
        (sym==:PV₊electrical_control₊Q_lim) && continue
        (sym==:PV₊electrical_control₊ΔQ) && continue
        (sym==:PV₊electrical_control₊s_Q) && continue
        (sym==:PV₊electrical_control₊s_Qint) && continue
        (sym==:PV₊electrical_control₊V_in) && continue
        (sym==:PV₊electrical_control₊V_lima) && continue
        (sym==:PV₊electrical_control₊V_con) && continue
        (sym==:PV₊electrical_control₊V_limb) && continue
        (sym==:PV₊electrical_control₊ΔV) && continue
        (sym==:PV₊electrical_control₊s_V) && continue
        (sym==:PV₊electrical_control₊s_Vint) && continue
        (sym==:PV₊electrical_control₊I_in) && continue
        (sym==:PV₊electrical_control₊I_lim) && continue
        (sym==:PV₊electrical_control₊I_t) && continue
        (sym==:PV₊electrical_control₊ΔI) && continue
        (sym==:PV₊electrical_control₊I_qin) && continue
        (sym==:PV₊electrical_control₊I_qcon) && continue
        (sym==:PV₊electrical_control₊I_sum) && continue
        (sym==:PV₊electrical_control₊I_qcmd) && continue
        (sym==:PV₊electrical_control₊P_refout) && continue
        (sym==:PV₊electrical_control₊P_lim) && continue
        (sym==:PV₊electrical_control₊ΔP) && continue
        (sym==:PV₊electrical_control₊ΔP_lim) && continue
        (sym==:PV₊electrical_control₊I_pref) && continue
        (sym==:PV₊electrical_control₊I_pcmd) && continue
        (sym==:PV₊electrical_control₊I_qmin) && continue
        (sym==:PV₊electrical_control₊I_qmax) && continue
        (sym==:PV₊electrical_control₊I_pmax) && continue
        (sym==:PV₊electrical_control₊I_pmin) && continue
        (sym==:PV₊converter_interface₊I_qrsum) && continue
        (sym==:PV₊converter_interface₊I_qrlim) && continue
        (sym==:PV₊converter_interface₊I_qr) && continue
        (sym==:PV₊converter_interface₊ΔV) && continue
        (sym==:PV₊converter_interface₊I_hv) && continue
        (sym==:PV₊converter_interface₊I_hvlim) && continue
        (sym==:PV₊converter_interface₊I_q) && continue
        (sym==:PV₊converter_interface₊ΔI_q) && continue
        (sym==:PV₊converter_interface₊ΔI_pr) && continue
        (sym==:PV₊converter_interface₊I_pr) && continue
        (sym==:PV₊converter_interface₊ΔI_prlim) && continue
        (sym==:PV₊converter_interface₊I_pg) && continue
        (sym==:PV₊converter_interface₊y) && continue
        (sym==:PV₊converter_interface₊I_p) && continue
        #(sym==:PV₊converter_interface₊V) && continue
        (sym==:PV₊converter_interface₊I_lvpl) && continue
        (sym==:PV₊V_t) && continue
        #(sym==:PV₊δ_v) && continue
        (sym==:PV₊pir) && continue
        (sym==:PV₊pii) && continue
        (sym==:PV₊pvr) && continue
        (sym==:PV₊pvi) && continue
        (sym==:PV₊P_gen) && continue
        (sym==:PV₊Q_gen) && continue
        (sym==:PV₊Vdiff) && continue
        (sym==:PV₊Vreg) && continue
        (sym==:PV₊Qbranch) && continue
        (sym==:PV₊Pbranch) && continue
        #= (sym==:BESS₊plant_control₊p_0) && continue
        (sym==:BESS₊plant_control₊Voltage_dip) && continue
        (sym==:BESS₊plant_control₊V_droop) && continue
        (sym==:BESS₊plant_control₊V_in) && continue
        (sym==:BESS₊plant_control₊V_fltr) && continue
        (sym==:BESS₊plant_control₊ΔV) && continue
        (sym==:BESS₊plant_control₊Q_fltr) && continue
        (sym==:BESS₊plant_control₊ΔQ) && continue
        (sym==:BESS₊plant_control₊ΔQ_in) && continue
        (sym==:BESS₊plant_control₊ΔQ_dbd) && continue
        (sym==:BESS₊plant_control₊Q_e) && continue
        (sym==:BESS₊plant_control₊Q_x) && continue
        (sym==:BESS₊plant_control₊Q_res) && continue
        (sym==:BESS₊plant_control₊Q_I) && continue
        (sym==:BESS₊plant_control₊Q_lim) && continue
        (sym==:BESS₊plant_control₊Q_ext) && continue
        (sym==:BESS₊plant_control₊Δf_deadband) && continue
        (sym==:BESS₊plant_control₊Δf_corr) && continue
        (sym==:BESS₊plant_control₊P_branchp) && continue
        (sym==:BESS₊plant_control₊f_e) && continue
        (sym==:BESS₊plant_control₊P_e) && continue
        (sym==:BESS₊plant_control₊P_lim) && continue
        (sym==:BESS₊plant_control₊P_refa) && continue
        (sym==:BESS₊plant_control₊P_ref) && continue
        (sym==:BESS₊electrical_control₊Voltage_dip) && continue
        (sym==:BESS₊electrical_control₊V_tfilt) && continue
        (sym==:BESS₊electrical_control₊V_tfiltlim) && continue
        (sym==:BESS₊electrical_control₊ΔV_t) && continue
        (sym==:BESS₊electrical_control₊ΔV_tdbd) && continue
        (sym==:BESS₊electrical_control₊I_qinj) && continue
        (sym==:BESS₊electrical_control₊P_PF) && continue
        (sym==:BESS₊electrical_control₊Q_con) && continue
        (sym==:BESS₊electrical_control₊Q_lim) && continue
        (sym==:BESS₊electrical_control₊ΔQ) && continue
        (sym==:BESS₊electrical_control₊s_Q) && continue
        (sym==:BESS₊electrical_control₊s_Qint) && continue
        (sym==:BESS₊electrical_control₊V_in) && continue
        (sym==:BESS₊electrical_control₊V_lima) && continue
        (sym==:BESS₊electrical_control₊V_con) && continue
        (sym==:BESS₊electrical_control₊V_limb) && continue
        (sym==:BESS₊electrical_control₊ΔV) && continue
        (sym==:BESS₊electrical_control₊s_V) && continue
        (sym==:BESS₊electrical_control₊s_Vint) && continue
        (sym==:BESS₊electrical_control₊I_in) && continue
        (sym==:BESS₊electrical_control₊I_lim) && continue
        (sym==:BESS₊electrical_control₊I_t) && continue
        (sym==:BESS₊electrical_control₊ΔI) && continue
        (sym==:BESS₊electrical_control₊I_qin) && continue
        (sym==:BESS₊electrical_control₊I_qcon) && continue
        (sym==:BESS₊electrical_control₊I_sum) && continue
        (sym==:BESS₊electrical_control₊I_qcmd) && continue
        (sym==:BESS₊electrical_control₊P_refout) && continue
        (sym==:BESS₊electrical_control₊P_lim) && continue
        (sym==:BESS₊electrical_control₊ΔP) && continue
        (sym==:BESS₊electrical_control₊ΔP_lim) && continue
        (sym==:BESS₊electrical_control₊I_pref) && continue
        (sym==:BESS₊electrical_control₊ΔI_p) && continue
        (sym==:BESS₊electrical_control₊I_pcmd) && continue
        (sym==:BESS₊electrical_control₊I_qmin) && continue
        (sym==:BESS₊electrical_control₊I_qmax) && continue
        (sym==:BESS₊electrical_control₊I_pmax) && continue
        (sym==:BESS₊electrical_control₊I_pmin) && continue
        (sym==:BESS₊electrical_control₊I_pmin_soc) && continue
        (sym==:BESS₊electrical_control₊I_pmax_soc) && continue
        (sym==:BESS₊electrical_control₊soc_Imin) && continue
        (sym==:BESS₊electrical_control₊soc_Imax) && continue
        (sym==:BESS₊electrical_control₊P_stor) && continue
        (sym==:BESS₊electrical_control₊soc) && continue
        (sym==:BESS₊electrical_control₊soc_lim) && continue
        (sym==:BESS₊electrical_control₊VDL1_out) && continue
        (sym==:BESS₊electrical_control₊VDL2_out) && continue
        (sym==:BESS₊converter_interface₊I_qrsum) && continue
        (sym==:BESS₊converter_interface₊I_qrlim) && continue
        (sym==:BESS₊converter_interface₊I_qr) && continue
        (sym==:BESS₊converter_interface₊ΔV) && continue
        (sym==:BESS₊converter_interface₊I_hv) && continue
        (sym==:BESS₊converter_interface₊I_hvlim) && continue
        (sym==:BESS₊converter_interface₊I_q) && continue
        (sym==:BESS₊converter_interface₊ΔI_q) && continue
        (sym==:BESS₊converter_interface₊ΔI_pr) && continue
        (sym==:BESS₊converter_interface₊I_pr) && continue
        (sym==:BESS₊converter_interface₊ΔI_prlim) && continue
        (sym==:BESS₊converter_interface₊I_pg) && continue
        (sym==:BESS₊converter_interface₊y) && continue
        (sym==:BESS₊converter_interface₊I_p) && continue
        (sym==:BESS₊converter_interface₊V) && continue
        (sym==:BESS₊converter_interface₊I_lvpl) && continue
        (sym==:BESS₊V_t) && continue
        (sym==:BESS₊δ_v) && continue
        (sym==:BESS₊pir) && continue
        (sym==:BESS₊pii) && continue
        (sym==:BESS₊pvr) && continue
        (sym==:BESS₊pvi) && continue
        (sym==:BESS₊PELEC) && continue
        (sym==:BESS₊P_gen) && continue
        (sym==:BESS₊Q_gen) && continue
        (sym==:BESS₊Vdiff) && continue
        (sym==:BESS₊Vreg) && continue
        (sym==:BESS₊Qbranch) && continue
        (sym==:BESS₊Pbranch) && continue =#
        #= (sym==:WT₊plant_control₊p_0) && continue
        (sym==:WT₊plant_control₊Voltage_dip) && continue
        (sym==:WT₊plant_control₊V_droop) && continue
        (sym==:WT₊plant_control₊V_in) && continue
        (sym==:WT₊plant_control₊V_fltr) && continue
        (sym==:WT₊plant_control₊ΔV) && continue
        (sym==:WT₊plant_control₊Q_fltr) && continue
        (sym==:WT₊plant_control₊ΔQ) && continue
        (sym==:WT₊plant_control₊ΔQ_in) && continue
        (sym==:WT₊plant_control₊ΔQ_dbd) && continue
        (sym==:WT₊plant_control₊Q_e) && continue
        (sym==:WT₊plant_control₊Q_x) && continue
        (sym==:WT₊plant_control₊Q_res) && continue
        (sym==:WT₊plant_control₊Q_I) && continue
        (sym==:WT₊plant_control₊Q_lim) && continue
        (sym==:WT₊plant_control₊Q_ext) && continue
        (sym==:WT₊plant_control₊Δf_deadband) && continue
        (sym==:WT₊plant_control₊Δf_corr) && continue
        (sym==:WT₊plant_control₊P_branchp) && continue
        (sym==:WT₊plant_control₊f_e) && continue
        (sym==:WT₊plant_control₊P_e) && continue
        (sym==:WT₊plant_control₊P_lim) && continue
        (sym==:WT₊plant_control₊P_refa) && continue
        (sym==:WT₊plant_control₊P_ref) && continue
        (sym==:WT₊electrical_control₊Voltage_dip) && continue
        (sym==:WT₊electrical_control₊V_tfilt) && continue
        (sym==:WT₊electrical_control₊V_tfiltlim) && continue
        (sym==:WT₊electrical_control₊ΔV_t) && continue
        (sym==:WT₊electrical_control₊ΔV_tdbd) && continue
        (sym==:WT₊electrical_control₊I_qinj) && continue
        (sym==:WT₊electrical_control₊P_PF) && continue
        (sym==:WT₊electrical_control₊Q_con) && continue
        (sym==:WT₊electrical_control₊Q_lim) && continue
        (sym==:WT₊electrical_control₊ΔQ) && continue
        (sym==:WT₊electrical_control₊s_Q) && continue
        (sym==:WT₊electrical_control₊V_in) && continue
        (sym==:WT₊electrical_control₊V_lima) && continue
        (sym==:WT₊electrical_control₊V_mod) && continue
        (sym==:WT₊electrical_control₊V_con) && continue
        (sym==:WT₊electrical_control₊V_limb) && continue
        (sym==:WT₊electrical_control₊ΔV) && continue
        (sym==:WT₊electrical_control₊s_V) && continue
        (sym==:WT₊electrical_control₊I_in) && continue
        (sym==:WT₊electrical_control₊I_lim) && continue
        (sym==:WT₊electrical_control₊I_t) && continue
        (sym==:WT₊electrical_control₊ΔI) && continue
        (sym==:WT₊electrical_control₊I_qin) && continue
        (sym==:WT₊electrical_control₊I_qcon) && continue
        (sym==:WT₊electrical_control₊I_sum) && continue
        (sym==:WT₊electrical_control₊I_qcmd) && continue
        (sym==:WT₊electrical_control₊P_in) && continue
        (sym==:WT₊electrical_control₊P_refout) && continue
        (sym==:WT₊electrical_control₊P_lim) && continue
        (sym==:WT₊electrical_control₊ΔP) && continue
        (sym==:WT₊electrical_control₊ΔP_lim) && continue
        (sym==:WT₊electrical_control₊I_pref) && continue
        (sym==:WT₊electrical_control₊I_pcmd) && continue
        (sym==:WT₊electrical_control₊I_qmin) && continue
        (sym==:WT₊electrical_control₊I_qmax) && continue
        (sym==:WT₊electrical_control₊I_pmax) && continue
        (sym==:WT₊electrical_control₊I_pmin) && continue
        (sym==:WT₊electrical_control₊I_pre) && continue
        (sym==:WT₊electrical_control₊I_post) && continue
        (sym==:WT₊electrical_control₊VDL1_out) && continue
        (sym==:WT₊electrical_control₊VDL2_out) && continue
        (sym==:WT₊converter_interface₊I_qrsum) && continue
        (sym==:WT₊converter_interface₊I_qrlim) && continue
        (sym==:WT₊converter_interface₊I_qr) && continue
        (sym==:WT₊converter_interface₊ΔV) && continue
        (sym==:WT₊converter_interface₊I_hv) && continue
        (sym==:WT₊converter_interface₊I_hvlim) && continue
        (sym==:WT₊converter_interface₊I_q) && continue
        (sym==:WT₊converter_interface₊ΔI_q) && continue
        (sym==:WT₊converter_interface₊ΔI_pr) && continue
        (sym==:WT₊converter_interface₊I_pr) && continue
        (sym==:WT₊converter_interface₊ΔI_prlim) && continue
        (sym==:WT₊converter_interface₊I_pg) && continue
        (sym==:WT₊converter_interface₊y) && continue
        (sym==:WT₊converter_interface₊I_p) && continue
        (sym==:WT₊converter_interface₊V) && continue
        (sym==:WT₊converter_interface₊I_lvpl) && continue
        (sym==:WT₊V_t) && continue
        (sym==:WT₊δ_v) && continue
        (sym==:WT₊pir) && continue
        (sym==:WT₊pii) && continue
        (sym==:WT₊WTr) && continue
        (sym==:WT₊WTi) && continue
        (sym==:WT₊P_gen) && continue
        (sym==:WT₊Q_gen) && continue
        (sym==:WT₊Vdiff) && continue
        (sym==:WT₊Vreg) && continue
        (sym==:WT₊Qbranch) && continue
        (sym==:WT₊Pbranch) && continue
        (sym==:WT₊drive_train₊w_tin) && continue
        (sym==:WT₊drive_train₊w_add) && continue
        (sym==:WT₊drive_train₊w_tsum) && continue
        (sym==:WT₊drive_train₊w_t) && continue
        (sym==:WT₊drive_train₊w_gin) && continue
        (sym==:WT₊drive_train₊Δw_g) && continue
        (sym==:WT₊drive_train₊w_g) && continue
        (sym==:WT₊drive_train₊w_gint) && continue
        (sym==:WT₊drive_train₊Δw) && continue =#
        set_default!(bus1, sym, get_guess(bus1, sym))
    end

    s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol, nwtol)
    #dump_initial_state(bus1)
    init_residual(bus1; verbose=true)

    prob = ODEProblem(nw, uflat(s0), (0,5), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end

# Default GENROE machine for SMIB system for testing controllers
function default_controller_smib_genroe()
    # Machine parameters from OpenIPSL IEEET1 test case (lines 5-25)
    S_b = 100e6
    M_b = 100e6
    H = 4.28
    D = 0
    # V_b = 400e3
    # ω_b = 2π*50

    # GENROE machine parameters (matching OpenIPSL test exactly)
    Tpd0 = 5
    Tppd0 = 0.07
    Tpq0 = 0.9
    Tppq0 = 0.09
    Xd = 1.84
    Xq = 1.75
    Xpd = 0.41
    Xpq = 0.6
    Xppd = 0.2
    Xppq = 0.2
    Xl = 0.12
    S10 = 0.11
    S12 = 0.39
    R_a = 0
    # angle_0 = 0.070492225331847
    # P_0 = 40000000
    # Q_0 = 5416582
    # v_0 = 1

    PSSE_GENROE(;
        Tpd0, Tppd0, Tpq0, Tppq0, H, D,
        Xd, Xq, Xpd, Xpq, Xppd, Xppq, Xl,
        S10, S12, R_a,
        M_b, S_b,
        pmech_input=false,
        efd_input=true,
        name=:machine
    )
end

function ref_rms_error(sol, csv, idx, col)
    t = csv[:, "time"]
    _ref = csv[:, col]

    # we find alle t in ti, where the jump occurs
    # at these points, we ignore the error because its hard to define with left/right limit
    Δtmean = mean(diff(t))
    _ti_jump = findall(Δ -> Δ < Δtmean/10000, diff(t))
    ti_jump = unique!(sort!(vcat(_ti_jump, _ti_jump .+ 1)))

    deleteat!(t, ti_jump)
    deleteat!(_ref, ti_jump)

    _sim = sol(t, idxs=idx).u

    norm(_ref .- _sim) / sqrt(length(_ref))
end