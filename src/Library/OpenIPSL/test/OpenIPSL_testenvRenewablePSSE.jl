function OpenIPSL_RePSSE_pv(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
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
        (sym==:PV₊repca₊p_0) && continue
        (sym==:PV₊repca₊Voltage_dip) && continue
        (sym==:PV₊repca₊V_droop) && continue
        (sym==:PV₊repca₊V_in) && continue
        (sym==:PV₊repca₊V_fltr) && continue
        (sym==:PV₊repca₊ΔV) && continue
        (sym==:PV₊repca₊Q_fltr) && continue
        (sym==:PV₊repca₊ΔQ) && continue
        (sym==:PV₊repca₊ΔQ_in) && continue
        (sym==:PV₊repca₊ΔQ_dbd) && continue
        (sym==:PV₊repca₊Q_e) && continue
        (sym==:PV₊repca₊Q_x) && continue
        (sym==:PV₊repca₊Q_res) && continue
        (sym==:PV₊repca₊Q_I) && continue
        (sym==:PV₊repca₊Q_lim) && continue
        (sym==:PV₊repca₊Q_ext) && continue
        (sym==:PV₊repca₊Δf_deadband) && continue
        (sym==:PV₊repca₊Δf_corr) && continue
        (sym==:PV₊repca₊P_branchp) && continue
        (sym==:PV₊repca₊f_e) && continue
        (sym==:PV₊repca₊P_e) && continue
        (sym==:PV₊repca₊P_lim) && continue
        (sym==:PV₊repca₊P_refa) && continue
        (sym==:PV₊repca₊P_ref) && continue
        (sym==:PV₊reecb₊Voltage_dip) && continue
        (sym==:PV₊reecb₊V_tfilt) && continue
        (sym==:PV₊reecb₊V_tfiltlim) && continue
        (sym==:PV₊reecb₊ΔV_t) && continue
        (sym==:PV₊reecb₊ΔV_tdbd) && continue
        (sym==:PV₊reecb₊I_qinj) && continue
        (sym==:PV₊reecb₊P_PF) && continue
        (sym==:PV₊reecb₊Q_con) && continue
        (sym==:PV₊reecb₊Q_lim) && continue
        (sym==:PV₊reecb₊ΔQ) && continue
        (sym==:PV₊reecb₊s_Q) && continue
        (sym==:PV₊reecb₊s_Qint) && continue
        (sym==:PV₊reecb₊V_in) && continue
        (sym==:PV₊reecb₊V_lima) && continue
        (sym==:PV₊reecb₊V_con) && continue
        (sym==:PV₊reecb₊V_limb) && continue
        (sym==:PV₊reecb₊ΔV) && continue
        (sym==:PV₊reecb₊s_V) && continue
        (sym==:PV₊reecb₊s_Vint) && continue
        (sym==:PV₊reecb₊I_in) && continue
        (sym==:PV₊reecb₊I_lim) && continue
        (sym==:PV₊reecb₊I_t) && continue
        (sym==:PV₊reecb₊ΔI) && continue
        (sym==:PV₊reecb₊I_qin) && continue
        (sym==:PV₊reecb₊I_qcon) && continue
        (sym==:PV₊reecb₊I_sum) && continue
        (sym==:PV₊reecb₊I_qcmd) && continue
        (sym==:PV₊reecb₊P_refout) && continue
        (sym==:PV₊reecb₊P_lim) && continue
        (sym==:PV₊reecb₊ΔP) && continue
        (sym==:PV₊reecb₊ΔP_lim) && continue
        (sym==:PV₊reecb₊I_pref) && continue
        (sym==:PV₊reecb₊I_pcmd) && continue
        (sym==:PV₊reecb₊I_qmin) && continue
        (sym==:PV₊reecb₊I_qmax) && continue
        (sym==:PV₊reecb₊I_pmax) && continue
        (sym==:PV₊reecb₊I_pmin) && continue
        (sym==:PV₊regca₊I_qrsum) && continue
        (sym==:PV₊regca₊I_qrlim) && continue
        (sym==:PV₊regca₊I_qr) && continue
        (sym==:PV₊regca₊ΔV) && continue
        (sym==:PV₊regca₊I_hv) && continue
        (sym==:PV₊regca₊I_hvlim) && continue
        (sym==:PV₊regca₊I_q) && continue
        (sym==:PV₊regca₊ΔI_q) && continue
        (sym==:PV₊regca₊ΔI_pr) && continue
        (sym==:PV₊regca₊I_pr) && continue
        (sym==:PV₊regca₊ΔI_prlim) && continue
        (sym==:PV₊regca₊I_pg) && continue
        (sym==:PV₊regca₊y) && continue
        (sym==:PV₊regca₊I_p) && continue
        #(sym==:PV₊regca₊V) && continue
        (sym==:PV₊regca₊I_lvpl) && continue
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


function OpenIPSL_RePSSE_bess(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
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

    #=
    for sym in sym(bus1)
        has_guess(bus1, sym) || continue
        (sym==:BESS₊repca₊p_0) && continue
        (sym==:BESS₊repca₊Voltage_dip) && continue
        (sym==:BESS₊repca₊V_droop) && continue
        (sym==:BESS₊repca₊V_in) && continue
        (sym==:BESS₊repca₊V_fltr) && continue
        (sym==:BESS₊repca₊ΔV) && continue
        (sym==:BESS₊repca₊Q_fltr) && continue
        (sym==:BESS₊repca₊ΔQ) && continue
        (sym==:BESS₊repca₊ΔQ_in) && continue
        (sym==:BESS₊repca₊ΔQ_dbd) && continue
        (sym==:BESS₊repca₊Q_e) && continue
        (sym==:BESS₊repca₊Q_x) && continue
        (sym==:BESS₊repca₊Q_res) && continue
        (sym==:BESS₊repca₊Q_I) && continue
        (sym==:BESS₊repca₊Q_lim) && continue
        (sym==:BESS₊repca₊Q_ext) && continue
        (sym==:BESS₊repca₊Δf_deadband) && continue
        (sym==:BESS₊repca₊Δf_corr) && continue
        (sym==:BESS₊repca₊P_branchp) && continue
        (sym==:BESS₊repca₊f_e) && continue
        (sym==:BESS₊repca₊P_e) && continue
        (sym==:BESS₊repca₊P_lim) && continue
        (sym==:BESS₊repca₊P_refa) && continue
        (sym==:BESS₊repca₊P_ref) && continue
        (sym==:BESS₊reecc₊Voltage_dip) && continue
        (sym==:BESS₊reecc₊V_tfilt) && continue
        (sym==:BESS₊reecc₊V_tfiltlim) && continue
        (sym==:BESS₊reecc₊ΔV_t) && continue
        (sym==:BESS₊reecc₊ΔV_tdbd) && continue
        (sym==:BESS₊reecc₊I_qinj) && continue
        (sym==:BESS₊reecc₊P_PF) && continue
        (sym==:BESS₊reecc₊Q_con) && continue
        (sym==:BESS₊reecc₊Q_lim) && continue
        (sym==:BESS₊reecc₊ΔQ) && continue
        (sym==:BESS₊reecc₊s_Q) && continue
        (sym==:BESS₊reecc₊s_Qint) && continue
        (sym==:BESS₊reecc₊V_in) && continue
        (sym==:BESS₊reecc₊V_lima) && continue
        (sym==:BESS₊reecc₊V_con) && continue
        (sym==:BESS₊reecc₊V_limb) && continue
        (sym==:BESS₊reecc₊ΔV) && continue
        (sym==:BESS₊reecc₊s_V) && continue
        (sym==:BESS₊reecc₊s_Vint) && continue
        (sym==:BESS₊reecc₊I_in) && continue
        (sym==:BESS₊reecc₊I_lim) && continue
        (sym==:BESS₊reecc₊I_t) && continue
        (sym==:BESS₊reecc₊ΔI) && continue
        (sym==:BESS₊reecc₊I_qin) && continue
        (sym==:BESS₊reecc₊I_qcon) && continue
        (sym==:BESS₊reecc₊I_sum) && continue
        (sym==:BESS₊reecc₊I_qcmd) && continue
        (sym==:BESS₊reecc₊P_refout) && continue
        (sym==:BESS₊reecc₊P_lim) && continue
        (sym==:BESS₊reecc₊ΔP) && continue
        (sym==:BESS₊reecc₊ΔP_lim) && continue
        (sym==:BESS₊reecc₊I_pref) && continue
        (sym==:BESS₊reecc₊ΔI_p) && continue
        (sym==:BESS₊reecc₊I_pcmd) && continue
        (sym==:BESS₊reecc₊I_qmin) && continue
        (sym==:BESS₊reecc₊I_qmax) && continue
        (sym==:BESS₊reecc₊I_pmax) && continue
        (sym==:BESS₊reecc₊I_pmin) && continue
        (sym==:BESS₊reecc₊I_pmin_soc) && continue
        (sym==:BESS₊reecc₊I_pmax_soc) && continue
        (sym==:BESS₊reecc₊soc_Imin) && continue
        (sym==:BESS₊reecc₊soc_Imax) && continue
        (sym==:BESS₊reecc₊P_stor) && continue
        (sym==:BESS₊reecc₊soc) && continue
        (sym==:BESS₊reecc₊soc_lim) && continue
        (sym==:BESS₊reecc₊VDL1_out) && continue
        (sym==:BESS₊reecc₊VDL2_out) && continue
        (sym==:BESS₊regca₊I_qrsum) && continue
        (sym==:BESS₊regca₊I_qrlim) && continue
        (sym==:BESS₊regca₊I_qr) && continue
        (sym==:BESS₊regca₊ΔV) && continue
        (sym==:BESS₊regca₊I_hv) && continue
        (sym==:BESS₊regca₊I_hvlim) && continue
        (sym==:BESS₊regca₊I_q) && continue
        (sym==:BESS₊regca₊ΔI_q) && continue
        (sym==:BESS₊regca₊ΔI_pr) && continue
        (sym==:BESS₊regca₊I_pr) && continue
        (sym==:BESS₊regca₊ΔI_prlim) && continue
        (sym==:BESS₊regca₊I_pg) && continue
        (sym==:BESS₊regca₊y) && continue
        (sym==:BESS₊regca₊I_p) && continue
        (sym==:BESS₊regca₊V) && continue
        (sym==:BESS₊regca₊I_lvpl) && continue
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
        (sym==:BESS₊Pbranch) && continue
        set_default!(bus1, sym, get_guess(bus1, sym))
    end
=#

    s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol, nwtol)
    #dump_initial_state(bus1)
    init_residual(bus1; verbose=true)

    prob = ODEProblem(nw, uflat(s0), (0,5), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end

function OpenIPSL_RePSSE_wt(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
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
   #=
    for sym in sym(bus1)
        has_guess(bus1, sym) || continue
        #=(sym==:WT₊repca₊p_0) && continue
        (sym==:WT₊repca₊Voltage_dip) && continue
        (sym==:WT₊repca₊V_droop) && continue
        (sym==:WT₊repca₊V_in) && continue
        (sym==:WT₊repca₊V_fltr) && continue
        (sym==:WT₊repca₊ΔV) && continue
        (sym==:WT₊repca₊Q_fltr) && continue
        (sym==:WT₊repca₊ΔQ) && continue
        (sym==:WT₊repca₊ΔQ_in) && continue
        (sym==:WT₊repca₊ΔQ_dbd) && continue
        (sym==:WT₊repca₊Q_e) && continue
        (sym==:WT₊repca₊Q_x) && continue
        (sym==:WT₊repca₊Q_res) && continue
        (sym==:WT₊repca₊Q_I) && continue
        (sym==:WT₊repca₊Q_lim) && continue
        (sym==:WT₊repca₊Q_ext) && continue
        (sym==:WT₊repca₊Δf_deadband) && continue
        (sym==:WT₊repca₊Δf_corr) && continue
        (sym==:WT₊repca₊P_branchp) && continue
        (sym==:WT₊repca₊f_e) && continue
        (sym==:WT₊repca₊P_e) && continue
        (sym==:WT₊repca₊P_lim) && continue
        #(sym==:WT₊repca₊P_refa) && continue
        (sym==:WT₊repca₊P_ref) && continue
        (sym==:WT₊reeca₊Voltage_dip) && continue
        (sym==:WT₊reeca₊V_tfilt) && continue
        (sym==:WT₊reeca₊V_tfiltlim) && continue
        (sym==:WT₊reeca₊ΔV_t) && continue
        (sym==:WT₊reeca₊ΔV_tdbd) && continue
        (sym==:WT₊reeca₊I_qinj) && continue
        (sym==:WT₊reeca₊P_PF) && continue
        (sym==:WT₊reeca₊Q_con) && continue
        (sym==:WT₊reeca₊Q_lim) && continue
        (sym==:WT₊reeca₊ΔQ) && continue
        (sym==:WT₊reeca₊s_Q) && continue
        (sym==:WT₊reeca₊V_in) && continue
        (sym==:WT₊reeca₊V_lima) && continue
        (sym==:WT₊reeca₊V_mod) && continue
        (sym==:WT₊reeca₊V_con) && continue
        (sym==:WT₊reeca₊V_limb) && continue
        (sym==:WT₊reeca₊ΔV) && continue
        (sym==:WT₊reeca₊s_V) && continue
        (sym==:WT₊reeca₊I_in) && continue
        (sym==:WT₊reeca₊I_lim) && continue
        (sym==:WT₊reeca₊I_t) && continue
        (sym==:WT₊reeca₊ΔI) && continue
        (sym==:WT₊reeca₊I_qin) && continue
        (sym==:WT₊reeca₊I_qcon) && continue
        (sym==:WT₊reeca₊I_sum) && continue
        (sym==:WT₊reeca₊I_qcmd) && continue
        (sym==:WT₊reeca₊P_in) && continue
        #(sym==:WT₊reeca₊P_refout) && continue
        (sym==:WT₊reeca₊P_lim) && continue
        (sym==:WT₊reeca₊ΔP) && continue
        (sym==:WT₊reeca₊ΔP_lim) && continue
        (sym==:WT₊reeca₊I_pref) && continue
        (sym==:WT₊reeca₊I_pcmd) && continue
        (sym==:WT₊reeca₊I_qmin) && continue
        (sym==:WT₊reeca₊I_qmax) && continue
        (sym==:WT₊reeca₊I_pmax) && continue
        (sym==:WT₊reeca₊I_pmin) && continue
        (sym==:WT₊reeca₊I_pre) && continue
        (sym==:WT₊reeca₊I_post) && continue
        (sym==:WT₊reeca₊VDL1_out) && continue
        (sym==:WT₊reeca₊VDL2_out) && continue
        (sym==:WT₊regca₊I_qrsum) && continue
        (sym==:WT₊regca₊I_qrlim) && continue
        (sym==:WT₊regca₊I_qr) && continue
        (sym==:WT₊regca₊ΔV) && continue
        (sym==:WT₊regca₊I_hv) && continue
        (sym==:WT₊regca₊I_hvlim) && continue
        (sym==:WT₊regca₊I_q) && continue
        (sym==:WT₊regca₊ΔI_q) && continue
        (sym==:WT₊regca₊ΔI_pr) && continue
        (sym==:WT₊regca₊I_pr) && continue
        (sym==:WT₊regca₊ΔI_prlim) && continue
        #(sym==:WT₊regca₊I_pg) && continue
        (sym==:WT₊regca₊y) && continue
        (sym==:WT₊regca₊I_p) && continue
        #(sym==:WT₊regca₊V) && continue
        (sym==:WT₊regca₊I_lvpl) && continue
        (sym==:WT₊V_t) && continue
        #(sym==:WT₊δ_v) && continue
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
    =#

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


function OpenIPSL_RePSSE_pv_pf(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
    # copy constructor and set vidxs
    bus1 = VertexModel(_bus1, vidx=1, name=:GEN1)
    @named junction = compile_bus(MTKBus(), vidx=2)
    loopback = LoopbackConnection(; src=:GEN1, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])
    #v_0 = 1.0
    #angle_0 = 0

    #@named slack =  SlackDifferential()
    #busmodel = MTKBus(slack; name=:slack_src)
    #bus2 = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0), vidx=2)
    slack = compile_bus(SlackAlgebraic(name=:slack_src), vidx=3)

    # line
    U_b = 230000
    S_b = 100000000
    Z_b = U_b^2/S_b
    pwLine = MTKLine(PiLine(; name=:PwLine))
    line = compile_line(pwLine; name=:pwLine,
        src=:junction, dst=:slack_src,
        PwLine₊X=1/Z_b, PwLine₊R=1/Z_b)

    buses = [bus1, junction, slack]
    lines = [line, loopback]
    nw = Network(buses, lines; warn_order=false)

    verbose = true
    pfnw=nothing
    pfs0=nothing
    pfs=nothing
    pfnw = isnothing(pfnw) ? powerflow_model(nw) : pfnw
    pfs0 = isnothing(pfs0) ? NWState(pfnw) : pfnw
    pfs = solve_powerflow(nw; pfnw, pfs0, verbose,t=0)
    println(show_powerflow(pfs))
    interface_vals = interface_values(pfs)
    println(interface_vals)
    # pfnw = powerflow_model(nw)
    # pfs = solve_powerflow(pfnw)

    if just_init
        s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol=Inf, nwtol=Inf)
        return s0
    end

    s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol, nwtol)
    #dump_initial_state(bus1)
    init_residual(bus1; verbose=true)

    prob = ODEProblem(nw, uflat(s0), (0,5), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end
