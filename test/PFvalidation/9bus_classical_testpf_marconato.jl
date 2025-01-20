using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks
using CairoMakie
using CSV
using DataFrames

@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset, Vset=nothing)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end

@mtkmodel ClassicBus begin
    @components begin
        machine = Library.ClassicalMachine_powerfactory(;
            p_m_input=false,
            S_b=100,
            V_b=18,
            ω_b=2π*60,
            X′_d,
            R_s=0.0026,
            vf_set=nothing,
            p_m_set=nothing,
            H,
             D=0
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

#aus Doku:
#Ae, Be = Library.solve_ceilf(3.3 => 0.6602, 4.5 => 4.2662)
@mtkmodel SauerPaiBus begin
    @components begin
        machine = Library.SauerPaiMachine(;
            S_b=247.5,
            V_b=16.5,
            ω_b=2π*60,
            vf_input=false,
            τ_m_input=false,
            H,
            R_s=0, 
            X_d=0.36135, 
            X_q=0.2398275, 
            X′_d=0.15048, 
            X′_q=0, 
            X″_d=0.1,
            X″_q=0.1, 
            X_ls=0.08316,
            T′_d0=8.96,
            T″_d0=0.075,
            T′_q0=0, 
            T″_q0=0.15
        )
        #avr = AVRTypeI(vr_min=-5, vr_max=5, Ka=20, Ta=0.2, Kf=0.063,
            #Tf=0.35, Ke=1, Te=0.314, Ae, Be, tmeas_lag=false)
        #gov = TGOV1(R=0.05, T1=0.05, T2=2.1, T3=7.0, DT=0, V_max=5, V_min=-5)
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
        #connect(machine.v_mag_out, avr.vh)
        #connect(avr.vf, machine.vf_in)
        #connect(gov.τ_m, machine.τ_m_in)
        #connect(machine.ωout, gov.ω_meas)
    end
end

@mtkmodel MarconatoBus begin
    @components begin
        machine = Library.MarconatoModel(;
            S_b,
            V_b,
            Ω_b=2π*60,
            vf_input=false,
            τ_m_input=false,
            H,
            D=0,
            R_s, 
            X_d, 
            X_q, 
            X′_d, 
            X′_q, 
            X″_d,
            X″_q, 
            X_ls,
            T′_d0,
            T″_d0,
            T′_q0, 
            T″_q0,
            T_AA
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end 

# generate all MTK bus models
@named mtkbus1 = MarconatoBus(; machine__S_b=247.5, machine__V_b=16.5, machine__H=23.64, machine__R_s=0, machine__X_d=0.36135, machine__X_q=0.2398275, machine__X′_d=0.15048, machine__X′_q=0.15, machine__X″_d=0.1, machine__X″_q=0.1, machine__X_ls=0.08316, machine__T′_d0=8.96, machine__T″_d0=0.075, machine__T′_q0=0.535, machine__T″_q0=0.15, machine__T_AA=0) #T′_q0, X′_q gibts im eigentlichen Modell nicht in PF, da dort Schenkelpolläufer -> Werte von hier verwendet!
@named mtkbus2 = MarconatoBus(; machine__S_b=192, machine__V_b=18, machine__H=6.40, machine__R_s=0.005, machine__X_d=1.719936, machine__X_q=1.65984, machine__X′_d=0.230016, machine__X′_q=0.378084, machine__X″_d=0.2, machine__X″_q=0.2, machine__X_ls=0.100032, machine__T′_d0=6, machine__T″_d0=0.0575, machine__T′_q0=0.535, machine__T″_q0=0.0945, machine__T_AA=0)
@named mtkbus3 = MarconatoBus(; machine__S_b=128, machine__V_b=13.8, machine__H=3.01, machine__R_s=0, machine__X_d=1.68, machine__X_q=1.609984, machine__X′_d=0.232064, machine__X′_q=0.32, machine__X″_d=0.2, machine__X″_q=0.2, machine__X_ls=0.094976, machine__T′_d0=5.89, machine__T″_d0=0.0575, machine__T′_q0=0.6, machine__T″_q0=0.08, machine__T_AA=0)
#@named mtkbus2 = ClassicBus(; machine__H= 6.40, machine__X′_d=0.1198)
#@named mtkbus3 = ClassicBus(; machine__H= 3.01, machine__X′_d=0.1813)
@named mtkbus4 = MTKBus()
@named mtkbus5 = LoadBus(;load__Pset=-1.25, load__Qset=-0.5)
@named mtkbus6 = LoadBus(;load__Pset=-0.90, load__Qset=-0.3)
@named mtkbus7 = MTKBus()
@named mtkbus8 = LoadBus(;load__Pset=-1.0, load__Qset=-0.35)
@named mtkbus9 = MTKBus()


# generate the dynamic component functions
@named bus1 = Bus(mtkbus1; vidx=1, pf=pfSlack(V=1.04)) 
@named bus2 = Bus(mtkbus2; vidx=2, pf=pfPV(V=1.025, P=1.63))
@named bus3 = Bus(mtkbus3; vidx=3, pf=pfPV(V=1.025, P=0.85))
@named bus4 = Bus(mtkbus4; vidx=4)
@named bus5 = Bus(mtkbus5; vidx=5, pf=pfPQ(P=-1.25, Q=-0.5))
@named bus6 = Bus(mtkbus6; vidx=6, pf=pfPQ(P=-0.9, Q=-0.3))
@named bus7 = Bus(mtkbus7; vidx=7)
@named bus8 = Bus(mtkbus8; vidx=8, pf=pfPQ(P=-1.0, Q=-0.35))
@named bus9 = Bus(mtkbus9; vidx=9)

# Branches
function piline(; R, X, B)
    @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
    MTKLine(pibranch)
end

function piline_shortcircuit(; R, X, B, pos, G_fault=0, B_fault=0)
    #faultimp = if (G_fault + B_fault) == 0
       # 0
    #else
        #1
    #end
    @named pibranch = PiLine_fault(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, G_fault, B_fault, pos)#, faultimp)
    MTKLine(pibranch)
end

function transformer(; R, X)
    @named transformer = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
    MTKLine(transformer)
end


@named l45 = Line(piline(; R=0.0100, X=0.0850, B=0.1760), src=4, dst=5)
@named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6)
@named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
@named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
@named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
@named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
@named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
@named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)
@named l57 = Line(piline_shortcircuit(; R=0.0320, X=0.1610, B=0.3060, pos=0.99), src=5, dst=7) #S_b = 100 MVA, U_b = 230 kV; 2 Ω also G_fault=0.003781

# build network
vertexfs = [bus1, bus2, bus3, bus4, bus5, bus6, bus7, bus8, bus9];
edgefs = [l45, l46, l69, l78, l89, t14, t27, t39, l57];
nw = Network(vertexfs, edgefs)

# solve powerflow and initialize
OpPoDyn.solve_powerflow!(nw)
OpPoDyn.initialize!(nw)

# get state for actual calculation
u0 = NWState(nw)


# create faults
affect1! = (integrator) -> begin
    if integrator.t == 0.0
        @info "Short circuit on line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊shortcircuit] = 1
    else
        error("Should not be reached.")
    end
end
cb_shortcircuit = PresetTimeCallback([0.0], affect1!) 

affect2! = (integrator) -> begin
    if integrator.t == 0.05
        @info "Deactivate line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊active] = 0
    else
        error("Should not be reached.")
    end
end
cb_deactivate = PresetTimeCallback([0.05], affect2!)

cb_set = CallbackSet(cb_shortcircuit, cb_deactivate)
prob = ODEProblem(nw, uflat(u0), (0,2), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P());
nothing

break # stop execution of script here

#### Machine Angle
ref = CSV.read("RotorAngle_standardModelPF_G1adjust.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Läuferwinkel (Marconato Machine Model)")
ts = range(sol.t[begin],sol.t[end],length=1000)
δ2 = sol(ts, idxs=VIndex(2, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
δ3 = sol(ts, idxs=VIndex(3, :machine₊δ)).u - sol(ts, idxs=VIndex(1, :machine₊δ)).u
lines!(ax, ts, rad2deg.(δ2); label="Bus 2")
lines!(ax, ref."Zeitpunkt in s", ref."firel in deg", color=Cycled(1), linestyle=:dash, label="Bus 2 ref")
lines!(ax, ts, rad2deg.(δ3); label="Bus 3")
lines!(ax, ref."Zeitpunkt in s", ref."firel in deg_1", color=Cycled(2), linestyle=:dash, label="Bus 3 ref")
axislegend(ax; position=:lt)
fig


#### Voltage Magnitude
ref = CSV.read("Bus5-7_standardModelPF_G1adjust.csv", DataFrame; header=2, decimal=',') #_faultimp #ref = CSV.read("Bus5-7_voltage.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude (Marconato Machine Model)")
ts = range(sol.t[begin],sol.t[end],length=1000)
umag5 = sqrt.(sol(ts; idxs=VIndex(5, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(5, :busbar₊u_i)).^2)
umag7 = sqrt.(sol(ts; idxs=VIndex(7, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(7, :busbar₊u_i)).^2)
lines!(ax, ts, umag5.u; label="Bus5")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u._1", color=Cycled(1), linestyle=:dash, label="Bus 5 ref")
lines!(ax, ts, umag7.u; label="Bus7")
lines!(ax, ref."Zeitpunkt in s", ref."u1, Betrag in p.u.", color=Cycled(2), linestyle=:dash, label="Bus 7 ref")
axislegend(ax; position=:rb)
fig


#=
# Plotting the Solution
fig = Figure(size=(1000,2000));
ax = Axis(fig[1, 1]; title="Active power")
for i in [1,2,3,5,6,8]
    lines!(ax, sol; idxs=VIndex(i,:busbar₊P), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[2, 1]; title="Voltage magnitude")
for i in 1:9
    lines!(ax, sol; idxs=VIndex(i,:busbar₊u_mag), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[3, 1]; title="Voltag angel")
for i in 1:9
    lines!(ax, sol; idxs=VIndex(i,:busbar₊u_arg), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
ax = Axis(fig[4, 1]; title="Frequency")
for i in 1:3
    lines!(ax, sol; idxs=VIndex(i,:machine₊ω), label="Bus $i", color=Cycled(i))
end
axislegend(ax)
fig
=#
