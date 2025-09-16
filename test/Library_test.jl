using OpPoDyn
using OpPoDyn.Library
using PowerDynamics.Library
using PowerDynamicsTesting
PowerDynamicsTesting.set_reference_dir(OpPoDyn)
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: ModelingToolkit as MTK
using OrderedCollections
using DiffEqCallbacks
using CairoMakie
using Test

@info "Start Library tests"

@testset "DynawoPiLine" begin
    @named branchA = DynawoPiLine(XPu=0.022522)
    @named branchB = DynawoPiLine(XPu=0.04189)
    line = Line(MTKLine(branchA, branchB));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_1" toi

    @named branchA = DynawoPiLine(XPu=0.022522, RPu=0.01)
    line = Line(MTKLine(branchA));
    toi = line_between_slacks(line)
    # isinteractive() && plottoi(toi)
    @reftest "DynawoPiLine_2" toi
end

@testset "Dynawo Machine test" begin
    # line model
    @named branchA = DynawoPiLine(XPu=0.022522)
    @named branchB = DynawoPiLine(XPu=0.04189)
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
    # system does not initialize properly?
    # maybe this would need a proper find fixpoint?
    @test_broken SciMLBase.successful_retcode(sol.retcode)

    #=
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
    =#
end

@testset "IPSLPSAT" begin
    @mtkmodel GenBus begin
        @components begin
            machine = OpPoDyn.Library.IPSLPSATOrder4(;
                Sn=100,
                Vn=18,
                V_b=18,
                ra=0,
                xd=0.8958,
                xq=0.8645,
                x1d=0.1198,
                x1q=0.1969,
                T1d0=6,
                T1q0=0.5350,
                M=12.8,
                D=0,
                ω_b=2π*50,
                S_b=100)
            avr = AVRTypeI(
                vref_input=true,
                vr_min=-5,
                vr_max=5,
                Ka=20,
                Ta=0.2,
                Kf=0.063,
                Tf=0.35,
                Ke=1,
                Te=0.314,
                Tr=0.001,
                A=0.0039,
                B=1.555)
            pmech = Blocks.Constant(k=1.63)
            vref = Blocks.Constant(k=1.120103884682511)
            busbar = BusBar()
        end
        @equations begin
            connect(vref.output, avr.vref)
            connect(avr.vf, machine.vf)
            connect(machine.v_mag_out, avr.v_mag)
            connect(pmech.output, machine.pm)
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
    isinteractive() && plottoi(toi)
    @reftest "IPSLSAT_on_slack" toi
end
