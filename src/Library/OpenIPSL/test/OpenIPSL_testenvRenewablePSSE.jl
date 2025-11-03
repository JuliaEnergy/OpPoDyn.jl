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
            p[:pwFault₊B] = 2
            p[:pwFault₊G] = 2
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