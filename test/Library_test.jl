using OpPoDyn
using OpPoDyn.Library
using OpPoDynTesting
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: ModelingToolkit as MTK
using OrderedCollections
using DiffEqCallbacks
isinteractive() && using GLMakie
using Test

@testset "DynawoPiLine" begin
    @named branchA = DynawoPiLine(X_pu=0.022522)
    @named branchB = DynawoPiLine(X_pu=0.04189)
    line = Line(MTKLine(branchA, branchB));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_1" toi

    @named branchA = DynawoPiLine(X_pu=0.022522, R_pu=0.01)
    line = Line(MTKLine(branchA));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_2" toi
end

@testset "PiLineTest" begin
    X = 0.1
    R = 0.2
    G = 0.3
    B = 0.4
    @named refbranch = DynawoPiLine(X_pu=X, R_pu=R, G_pu=G, B_pu=B)
    @named branch = PiLine(X=X, R=R, G_src=G, G_dst=G, B_src=B, B_dst=B)

    line1 = Line(MTKLine(refbranch))
    line2 = Line(MTKLine(branch))
    p1 = get_default.(Ref(line1), psym(line1))
    p2 = get_default.(Ref(line2), psym(line2))

    out_src_1 = zeros(2)
    out_dst_1 = zeros(2)
    out_src_2 = zeros(2)
    out_dst_2 = zeros(2)
    usrc = rand(2)
    udst = rand(2)
    line1.g(out_src_1, out_dst_1, usrc, udst, p1, NaN)
    line2.g(out_src_2, out_dst_2, usrc, udst, p2, NaN)
    @test out_src_1 ≈ out_src_2
    @test out_dst_1 ≈ out_dst_2
end

@testset "PiLine_shortcircuit" begin
    R = 0.3
    X = 0.16
    B = 0.3
    @named branch = PiLine_fault(R=R, X=X, B_src=B/2, B_dst=B/2)
    line = Line(MTKLine(branch));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "PiLine_faultTestA" toi

    @named branch = PiLine_fault(R=R, X=X, B_src=B/2, B_dst=B/2, sc=1)
    line = Line(MTKLine(branch));
    toi = line_between_slacks(line)
    @reftest "PiLine_faultTestB" toi

    @named branch = PiLine_fault(R=R, X=X, B_src=B/2, B_dst=B/2, sc=1, active=0)
    line = Line(MTKLine(branch));
    toi = line_between_slacks(line)
    @reftest "PiLine_faultTestC" toi
end

@testset "Swing bus" begin
    @named swing = Swing(P_m=1, D=0.1, M=0.005, θ=0, ω=1, V_mag=1)
    bus = Bus(MTKBus(swing));
    toi = bus_on_slack(bus)
    # isinteractive() && plottoi(toi)
    @reftest "SwingBus_1" toi

    # swing bus with load
    @named swing = Swing(P_m=1, D=0.1, M=0.005, θ=0, ω=1, V_mag=1)
    @named pqload = PQLoad(P_set=-0.5, Q_set=-0.2)
    bm = MTKBus(swing, pqload)
    @test length(full_equations(simplify_mtkbus(bm))) == 2
    bus = Bus(bm)
    toi = bus_on_slack(bus)
    toi["active power"]["electric power of swing"] = VIndex(2,:swing₊P_e)
    toi["active power"]["electric power of load"] = VIndex(2,:pqload₊P)
    toi["reactive power"]["electric power of load"] = VIndex(2,:pqload₊Q)
    # isinteractive() && plottoi(toi)
    @reftest "Swing_and_load" toi
end

@testset "Dynawo Machine test" begin
    # line model
    @named branchA = DynawoPiLine(X_pu=0.022522)
    @named branchB = DynawoPiLine(X_pu=0.04189)
    linem = MTKLine(branchA, branchB)
    linef = Line(linem);

    # genbus model
    @mtkmodel GenBus begin
        @components begin
            machine = DynawoMachine()
            # machine = Swing(Pm_input=true)
            excitation = Blocks.Constant(k=2.4659)
            pmech = Blocks.Constant(k=0.903)
            ωRef = Blocks.Constant(k=1.0)
            trafo = DynawoFixedRatioTransformer()
            busbar = BusBar()
        end
        @equations begin
            connect(machine.efdPu, excitation.output)
            connect(machine.PmPu, pmech.output)
            # connect(machine.Pm, pmech.output)
            connect(machine.ωRefPu, ωRef.output)
            connect(machine.terminal, trafo.dst)
            connect(trafo.src, busbar.terminal)
        end
    end
    @named genbus = GenBus()
    # genbus = pin_parameters(genbus)
    genf = Bus(genbus; verbose=false)

    @named slack = SlackDifferential(u_init_r=0.90081)
    slackf = Bus(slack)

    g = path_graph(2)
    nw = Network(g, [slackf, genf], linef)
    u0 = NWState(nw)
    u0.v[2, :machine₊θ]        = 1.2107
    u0.v[2, :machine₊λ_fPu]    = 1.1458
    u0.v[2, :machine₊λ_DPu]    = 0.89243
    u0.v[2, :machine₊λ_Q1Pu]   = -0.60044
    u0.v[2, :machine₊λ_Q2Pu]   = -0.60044
    u0.v[2, :machine₊ωPu]      = 1
    u0.v[2, :machine₊idPu]     = -0.91975
    u0.v[2, :machine₊iqPu]     = -0.39262
    u0.v[2, :machine₊MqSat′Pu] = 1.5292
    u0.v[2, :machine₊MdSat′Pu] = 1.5792
    u0.v[2, :machine₊ifPu]     = 1.4855
    u0.v[2, :machine₊iDPu]     = 0
    u0.v[2, :machine₊iQ2Pu]    = 0
    u0.v[2, :machine₊iQ1Pu]    = 0

    affect = function(int)
        p = NWParameter(int)
        p.v[2, :pmech₊k] = 0.923
        auto_dt_reset!(int)
        save_parameters!(int)
    end
    cb = PresetTimeCallback(10, affect)
    prob = ODEProblem(nw, uflat(u0), (0, 30), copy(pflat(u0)); callback=cb)
    sol = solve(prob, Rodas5P())

    plotspec = OrderedDict(
        "active power" => OrderedDict(
            "injection from bus" => VIndex(2, :busbar₊P),
            "PGenPu" => VIndex(2, :machine₊PGenPu)),
        "reactive power" => OrderedDict(
            "injection from bus" => VIndex(2, :busbar₊Q),
            "QGenPu" => VIndex(2, :machine₊QGenPu)),
        "voltage angle" => OrderedDict(
            "angle at bus" => VIndex(2, :busbar₊u_arg),
            "angle of rotor" => VIndex(2, :machine₊θ),
            "angle at slack" => VIndex(1, :busbar₊u_arg)),
        "voltage magnitude" => OrderedDict(
            "magnitude at bus" => VIndex(2, :busbar₊u_mag),
            "stator voltage" => VIndex(2, :machine₊UStatorPu),
            "magnitude at slack" => VIndex(1, :busbar₊u_mag)),
        "frequency" => OrderedDict(
            "frequency at machine" => VIndex(2, :machine₊ωPu)))
    toi = TrajectoriesOfInterest(sol, plotspec)
    isinteractive() && plottoi(toi)
    @reftest "dynawomachine_on_slack_prefstep" toi
end

@testset "IPSLPSAT" begin
    @mtkmodel GenBus begin
        @components begin
            machine = OpPoDyn.Library.IPSLPSATOrder4(;
                S_n=100,
                V_n=18,
                V_b=18,
                R_a=0,
                X_d=0.8958,
                X_q=0.8645,
                X′_d=0.1198,
                X′_q=0.1969,
                T′_d0=6,
                T′_q0=0.5350,
                M=12.8,
                D=0,
                ω_b=2π*50,
                S_b=100)
            avr = AVRTypeI(
                V_ref_input=true,
                V_rmin=-5,
                V_rmax=5,
                K_a=20,
                T_a=0.2,
                K_f=0.063,
                T_f=0.35,
                K_e=1,
                T_e=0.314,
                T_r=0.001,
                A=0.0039,
                B=1.555)
            pmech = Blocks.Constant(k=1.63)
            vref = Blocks.Constant(k=1.120103884682511)
            busbar = BusBar()
        end
        @equations begin
            connect(vref.output, avr.V_ref)
            connect(avr.V_f, machine.V_f)
            connect(machine.V_mag_out, avr.V_mag)
            connect(pmech.output, machine.P_m)
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()
    bus = Bus(mtkbus);

    # obtained from steadystate for now
    set_default!(bus, :busbar₊u_r, 1.0189261518036425)
    set_default!(bus, :busbar₊u_i, 0.06828069999522467)
    set_default!(bus, :busbar₊i_r, -1.6299999998860033)
    set_default!(bus, :busbar₊i_i, 0.4518059633240033)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "IPSLPSAT" toi
end


@testset "SauerPai Generator" begin
    @mtkmodel GenBus begin
        @components begin
            machine = OpPoDyn.Library.SauerPaiMachine(;
                V_f_input=false,
                τ_m_input=false,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X_d=0.146, X′_d=0.0608, X″_d=0.06,
                X_q=0.1, X′_q=0.0969, X″_q=0.06,
                R_s=0.000124,
                X_ls=0.01460,
                T′_d0=8.96,
                T″_d0=0.01,
                T′_q0=0.31,
                T″_q0=0.01,
                H=23.64,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()
    bus = Bus(mtkbus)

    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "SauerPai" toi
end

@testset "SauerPai Generator with AVR and GOV" begin
    A, B = Library.solve_ceilf(3.3 => 0.6602, 4.5 => 4.2662)

    @mtkmodel GenBus begin
        @components begin
            machine = OpPoDyn.Library.SauerPaiMachine(;
                V_f_input=true,
                τ_m_input=true,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X_d=0.146, X′_d=0.0608, X″_d=0.06,
                X_q=0.1, X′_q=0.0969, X″_q=0.06,
                R_s=0.000124,
                X_ls=0.01460,
                T′_d0=8.96,
                T″_d0=0.01,
                T′_q0=0.31,
                T″_q0=0.01,
                H=23.64,
            )
            avr = AVRTypeI(
                V_rmin=-5,
                V_rmax=5,
                K_a=20,
                T_a=0.2,
                K_f=0.063,
                T_f=0.35,
                K_e=1,
                T_e=0.314,
                A, B,
                t_meas_lag=false)
            gov = TGOV1(
                R=0.05,
                T_1=0.05,
                T_2=2.1,
                T_3=7.0,
                D=0,
                V_max=5,
                V_min=-5)
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
            connect(machine.V_mag_out, avr.V_mag)
            connect(avr.V_f, machine.V_f_in)
            connect(gov.τ_m, machine.τ_m_in)
            connect(machine.ω_out, gov.ω_meas)
        end
    end
    @named mtkbus = GenBus()

    bus = Bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "SauerPai_AVRGOV" toi
end

@testset "test loads" begin
    @named load = PQLoad(P_set=-0.5, Q_set=-0.5)
    bus = Bus(MTKBus(load));
    toi = bus_on_slack(bus)
    #isinteractive() && plottoi(toi)
    @reftest "PQLoad" toi

    @named load = VoltageDependentLoad(P_set=-0.5, Q_set=-0.5, V_n=1, α_P=1, α_Q=1)
    bus = Bus(MTKBus(load));
    toi = bus_on_slack(bus)
    #isinteractive() && plottoi(toi)
    @reftest "VolatageDependentLoad" toi

    @named load = ConstantYLoad(P_set=-0.5, Q_set=-0.5, V_set=1)
    bus = Bus(MTKBus(load));
    toi = bus_on_slack(bus)
    #isinteractive() && plottoi(toi)
    @reftest "ConstantYLoad" toi
end

@testset "test zip load" begin
    # constant power
    @named load = ZIPLoad(P_set=-1, Q_set=-1, V_set=1,
                          K_pZ=0, K_pI=0, K_pC=1,
                          K_qZ=0, K_qI=0, K_qC=1)
    bus = Bus(MTKBus(load))
    toi = bus_on_slack(bus)
    #isinteractive() && plottoi(toi)
    @reftest "ZIPLoad_powerconst" toi

    # constant current
    @named load = ZIPLoad(P_set=-1, Q_set=-1, V_set=1,
                          K_pZ=0, K_pI=1, K_pC=0,
                          K_qZ=0, K_qI=1, K_qC=0)
    bus = Bus(MTKBus(load))
    toi = bus_on_slack(bus)
    toi["current"] = OrderedDict(
        "current magnitude at bus" => VIndex(2,:busbar₊i_mag),
        "real current" => VIndex(2,:load₊terminal₊i_r),
        "imaginary current" => VIndex(2,:load₊terminal₊i_i),
    )
    #isinteractive() && plottoi(toi)
    @reftest "ZIPLoad_currentconst" toi

    # constant Z
    @named load = ZIPLoad(P_set=-1, Q_set=-1, V_set=1,
                          K_pZ=1, K_pI=0, K_pC=0,
                          K_qZ=1, K_qI=0, K_qC=0)
    bus = Bus(MTKBus(load))
    toi = bus_on_slack(bus)
    toi["current"] = OrderedDict(
        "current magnitude at bus" => VIndex(2,:busbar₊i_mag),
        "real current" => VIndex(2,:load₊terminal₊i_r),
        "imaginary current" => VIndex(2,:load₊terminal₊i_i),
    )
    #isinteractive() && plottoi(toi) # neither i nor p is constant
    @reftest "ZIPLoad_Zconst" toi
end

@testset "Classical machine" begin
    @mtkmodel GenBus begin
        @components begin
            machine = Library.ClassicalMachine(;
                τ_m_input=false,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X′_d=0.0608,
                R_s=0.000124,
                H=23.64,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()

    simp = simplify_mtkbus(mtkbus)
    full_equations(simp)
    observed(simp)

    bus = Bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "ClassicalMachine" toi
end

@testset "Classical machine PowerFactory" begin
    @mtkmodel GenBus begin
        @components begin
            machine = Library.ClassicalMachine_powerfactory(;
                P_m_input=false,
                S_b=100,
                V_b=18,
                ω_b=2π*60,
                X′_d=0.0608,
                R_s=0.000124,
                H=23.64,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()

    simp = simplify_mtkbus(mtkbus)
    full_equations(simp)
    observed(simp)

    bus = Bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "ClassicalMachine_PowerFactory" toi
end

@testset "PF Standard Model machine - salient pole" begin
    @mtkmodel GenBus begin
        @components begin
            machine = Library.StandardModel_pf(;
                τ_m_input=false,
                 V_f_input = false,
                S_b=100,
                S_n=100,
                V_b=18,
                V_n=18,
                ω_b=2π*60,
                salientpole=1,
                H=9,
                D=0,
                R_s=0,
                X_rld=0,
                X_rlq=0,
                X_d=0.36,
                X_q=0.24,
                X′_d=0.15,
                X″_d=0.1,
                X″_q=0.1,
                X_ls=0.08,
                T′_d0=9,
                T″_d0=0.07,
                T″_q0=0.15,
                dpe=0,
                dkd=0,
                cosn=1,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()

    simp = simplify_mtkbus(mtkbus)
    full_equations(simp)
    observed(simp)

    bus = Bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "PFStandardModel_salientpole" toi
end

@testset "PF Standard Model machine - cylindrical rotor" begin
    @mtkmodel GenBus begin
        @components begin
            machine = Library.StandardModel_pf(;
                τ_m_input=false,
                 V_f_input = false,
                S_b=100,
                S_n=100,
                V_b=18,
                V_n=18,
                ω_b=2π*60,
                salientpole=0,
                H=9,
                D=0,
                R_s=0,
                X_rld=0,
                X_rlq=0,
                X_d=0.36,
                X_q=0.24,
                X′_q=0.4,
                X′_d=0.15,
                X″_d=0.1,
                X″_q=0.1,
                X_ls=0.08,
                T′_d0=9,
                T′_q0=0.5,
                T″_d0=0.07,
                T″_q0=0.15,
                dpe=0,
                dkd=0,
                cosn=1,
            )
            busbar = BusBar()
        end
        @equations begin
            connect(machine.terminal, busbar.terminal)
        end
    end
    @named mtkbus = GenBus()

    simp = simplify_mtkbus(mtkbus)
    full_equations(simp)
    observed(simp)

    bus = Bus(mtkbus)
    set_voltage!(bus; mag=1.017, arg=0.0295)
    set_current!(bus; P=0.716, Q=0.3025)
    initialize_component!(bus)

    toi = bus_on_slack(bus; tmax=600, toilength=10_000)
    #isinteractive() && plottoi(toi)
    @reftest "PFStandardModel_nonsalientpole" toi
end

@testset "AVR model" begin
    E_1 = 3.5461
    E_2 = 4.7281
    S_e1 = 0.08
    S_e2 = 0.26
    se_quad = x -> Library.quadratic_ceiling(x, E_1, E_2, S_e1, S_e2)
    A, B = Library.solve_ceilf(E_1=>S_e1, E_2=>S_e2)
    se_exp  = x -> A* exp(B*x)

    if isinteractive()
        let fig = Figure()
            ax = Axis(fig[1,1])
            xs = range(0, 6; length=100)
            lines!(ax, xs, se_quad.(xs); label="quad")
            lines!(ax, xs, se_exp.(xs); label="exp")
            axislegend(ax)
            scatter!(ax, [(E_1, S_e1), (E_2, S_e2)])
            fig
        end
    end
end

