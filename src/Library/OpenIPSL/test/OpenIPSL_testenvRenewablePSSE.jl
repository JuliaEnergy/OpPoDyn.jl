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

    P_load = 30e6   # W
    Q_load = 50e6   # var (inductive)
    S_b    = 100e6
    U_b    = 230000
    Z_b    = U_b^2 / S_b

    function make_junction(G, B)
        @named junc_load = ConstantYLoad(G=G, B=B, allow_zero_conductance=true)
        busmodel = MTKBus(junc_load; name=:junction)
        junc = compile_bus(busmodel, vidx=2)
        _disable = ComponentAffect([], [:junc_load₊G, :junc_load₊B]) do u, p, ctx
            @info "Load disconnected at junction bus at t = $(ctx.t)s"
            p[:junc_load₊G] = 0
            p[:junc_load₊B] = 0
        end
        set_callback!(junc, PresetTimeComponentCallback(1.0, _disable))
        junc
    end

    function make_network(junction)
        loopback = LoopbackConnection(; src=:GEN1, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])
        #v_0 = 1.0
        #angle_0 = 0

        #@named slack =  SlackDifferential()
        #busmodel = MTKBus(slack; name=:slack_src)
        #bus2 = compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0), vidx=2)
        slack = compile_bus(SlackAlgebraic(name=:slack_src), vidx=3)

        # line (100× longer than original: R=X=100/Z_b)
        pwLine = MTKLine(PiLine(; name=:PwLine))
        line = compile_line(pwLine; name=:pwLine,
            src=:junction, dst=:slack_src,
            PwLine₊X=100/Z_b, PwLine₊R=100/Z_b,
            PwLine₊B_src=0.0, PwLine₊B_dst=0.0,
            PwLine₊G_src=0.0, PwLine₊G_dst=0.0)

        Network([bus1, junction, slack], [line, loopback]; warn_order=false)
    end

    # Iterative self-consistent admittance: G = P/(V_lf^2 * S_b) where V_lf is the power flow voltage
    # Converges when the power flow voltage no longer changes (typically 2-3 iterations)
    G_load, B_load = P_load/S_b, -Q_load/S_b
    V_lf = 1.0
    for _ in 1:50
        pfs_tmp = solve_powerflow(make_network(make_junction(G_load, B_load)); verbose=false)
        V_new   = pfs_tmp[VIndex(2, :busbar₊u_mag)]
        abs(V_new - V_lf) < 1e-10 && break
        V_lf   = V_new
        G_load = P_load / (V_lf^2 * S_b)
        B_load = -Q_load / (V_lf^2 * S_b)
    end
    @info "Self-consistent V_lf = $(round(V_lf; digits=6)) pu → G=$(round(G_load;digits=6)), B=$(round(B_load;digits=6))"

    # Pass 2: build network with self-consistent admittance
    junction = make_junction(G_load, B_load)
    nw       = make_network(junction)

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

    s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol, nwtol)
    #dump_initial_state(bus1)
    init_residual(bus1; verbose=true)

    prob = ODEProblem(nw, uflat(s0), (0,5), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end


function OpenIPSL_RePSSE_pv_pf_3bus(_bus1; ω_b=2π*50, just_init=false, tol=1e0, nwtol=1e0)
    bus1 = VertexModel(_bus1, vidx=1, name=:GEN1)
    @named junction  = compile_bus(MTKBus(), vidx=2)
    slack            = compile_bus(SlackAlgebraic(name=:slack_src), vidx=3)
    @named fault_bus = compile_bus(MTKBus(), vidx=4)

    loopback = LoopbackConnection(; src=:GEN1, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])

    U_b = 230000
    S_b = 100000000
    Z_b = U_b^2/S_b

    # Existing line: junction → slack_src (same parameters as 2-bus case)
    line_junc_slack = compile_line(
        MTKLine(PiLine(; name=:PwLine));
        name=:line_junc_slack,
        src=:junction, dst=:slack_src,
        PwLine₊X=1/Z_b, PwLine₊R=1/Z_b)

    # New line: fault_bus → junction (PiLine_fault, can be short-circuited)
    line_fault_junc = compile_line(
        MTKLine(PiLine_fault(; name=:PwLineFault));
        name=:line_fault_junc,
        src=:fault_bus, dst=:junction,
        PwLineFault₊X=1/Z_b, PwLineFault₊R=1/Z_b,
        PwLineFault₊G_src=0, PwLineFault₊B_src=0,
        PwLineFault₊G_dst=0, PwLineFault₊B_dst=0,
        PwLineFault₊pos=0.5)

    # t=0.1s: three-phase short circuit at 50% of line fault_bus→junction
    _enable_short = ComponentAffect([], [:PwLineFault₊shortcircuit]) do u, p, ctx
        @info "Short circuit activated on line $(ctx.src)→$(ctx.dst) at t = $(ctx.t)s"
        p[:PwLineFault₊shortcircuit] = 1
    end
    shortcircuit_cb = PresetTimeComponentCallback(0.1, _enable_short)

    # t=0.2s: disconnect line (fault clearing)
    _disable_line = ComponentAffect([], [:PwLineFault₊active]) do u, p, ctx
        @info "Line $(ctx.src)→$(ctx.dst) disconnected at t = $(ctx.t)s"
        p[:PwLineFault₊active] = 0
    end
    deactivate_cb = PresetTimeComponentCallback(0.2, _disable_line)

    set_callback!(line_fault_junc, (shortcircuit_cb, deactivate_cb))

    # New line: fault_bus → slack_src (same parameters as existing line)
    line_fault_slack = compile_line(
        MTKLine(PiLine(; name=:PwLine2));
        name=:line_fault_slack,
        src=:fault_bus, dst=:slack_src,
        PwLine2₊X=1/Z_b, PwLine2₊R=1/Z_b)

    buses = [bus1, junction, slack, fault_bus]
    lines = [loopback, line_junc_slack, line_fault_junc, line_fault_slack]
    nw = Network(buses, lines; warn_order=false)

    pfnw = powerflow_model(nw)
    pfs0 = NWState(pfnw)
    pfs  = solve_powerflow(nw; pfnw, pfs0, verbose=true)
    println(show_powerflow(pfs))
    println(interface_values(pfs))

    if just_init
        s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol=Inf, nwtol=Inf)
        return s0
    end

    s0 = initialize_from_pf!(nw; subverbose=[VIndex(1)], tol, nwtol)
    init_residual(bus1; verbose=true)

    prob = ODEProblem(nw, uflat(s0), (0,5), copy(pflat(s0)), callback=get_callbacks(nw))
    sol = solve(prob, Rodas5P())
    @assert SciMLBase.successful_retcode(sol) "Simulation was not successful: retcode=$(sol.retcode)"
    sol
end
