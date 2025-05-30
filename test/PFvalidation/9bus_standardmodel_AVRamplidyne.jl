using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using NetworkDynamicsInspector
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks
using CairoMakie
using CSV
using DataFrames
using Test

@mtkmodel LoadBus begin
    @components begin
        busbar = BusBar()
        load = ConstantYLoad(Pset, Qset, Vset=nothing)
    end
    @equations begin
        connect(load.terminal, busbar.terminal)
    end
end

@mtkmodel StandardBus begin
    @components begin
        machine = Library.StandardModel_pf(;
            S_b,
            Sn,
            V_b,
            Vn,
            ω_b=2π*60,
            H,
            D,
            vf_input=false,
            τ_m_input=false,
            R_s,
            X_rld,
            X_rlq,
            X_d,
            X_q,
            X′_q,
            X′_d,
            X″_d,
            X″_q,
            X_ls,
            T′_d0,
            T′_q0,
            T″_d0,
            T″_q0,
            dpe,
            dkd,
            cosn,
            salientpole,
            dpu,
            addmt,
            xmdm
        )
        busbar = BusBar()
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
    end
end

@mtkmodel StandardBus_AVR begin
    @components begin
        machine = Library.StandardModel_pf(;
            S_b,
            Sn,
            V_b,
            Vn,
            ω_b=2π*60,
            H,
            D,
            vf_input=true,
            τ_m_input=false,
            R_s,
            X_rld,
            X_rlq,
            X_d,
            X_q,
            X′_q,
            X′_d,
            X″_d,
            X″_q,
            X_ls,
            T′_d0,
            T″_d0,
            T′_q0,
            T″_q0,
            dpe,
            dkd,
            cosn,
            salientpole,
            dpu,
            addmt,
            xmdm
        )
        busbar = BusBar()
        avr = AVRTypeI(
            Ka, Ke, Kf,
            Ta, Tf, Te, Tr,
            vr_min, vr_max, anti_windup,
            A, B
        )
    end
    @equations begin
        connect(machine.terminal, busbar.terminal)
        connect(machine.v_mag_out, avr.v_mag)
        connect(avr.vf, machine.vf_in)
    end
end

# Branches
function piline(; R, X, B)
    @named pibranch = PiLine(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0)
    MTKLine(pibranch)
end

function piline_shortcircuit(; R, X, B, pos, G_fault=0, B_fault=0, faultimp=0)
    @named pibranch = PiLine_fault(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, G_fault, B_fault, pos, faultimp)
    MTKLine(pibranch)
end

function transformer(; R, X)
    @named transformer = PiLine(;R, X, B_src=0, B_dst=0, G_src=0, G_dst=0)
    MTKLine(transformer)
end

function lineparams_pu(; S_b, V_b, line_length, R_l_perkm, X_l_perkm, B_l_perkm, R_f_Ω=1, X_f_Ω=1)
    R_pu = (R_l_perkm * line_length) * S_b / V_b^2
    X_pu = (X_l_perkm * line_length) * S_b / V_b^2
    B_pu = (B_l_perkm * line_length) * V_b^2 / S_b #ω_b * 0.0095491* 10^(-6) #ω*C
    G_f_pu = 1/R_f_Ω * V_b^2 / S_b
    B_f_pu = 1/X_f_Ω * V_b^2 / S_b
    return (R=R_pu, X=X_pu, B=B_pu, G_f=G_f_pu, B_f=B_f_pu)
end

function trafoparams_pu(; S_b, S_n, X_pu_PF, R_pu_PF)
    X_pu_baseS = X_pu_PF * S_b/S_n
    R_pu_baseS = R_pu_PF * S_b/S_n
    return(R=R_pu_baseS, X=X_pu_baseS)
end



# generate all MTK bus models
primary_parameters_gen1 = Dict(
    Symbol("machine__", k) => v for (k, v) in
    (
        :S_b => 100e6,
        :Sn => 247.5e6,
        :V_b => 16.5,
        :Vn => 16.5,
        :ω_b => 2π*60,
        :H => 9.551516,
        :D => 0,
        :R_s => 0,
        :X_rld => 0,
        :X_rlq => 0,
        :X_d => 0.36135,
        :X_q => 0.2398275,
        :X′_d => 0.15048,
        :X′_q => 0.0001,
        :X″_d => 0.1,
        :X″_q => 0.1,
        :X_ls => 0.08316,
        :T′_d0 => 8.96,
        :T″_d0 => 0.075,
        :T′_q0 => 0.0001,
        :T″_q0 => 0.15,
        :cosn => 1,
        :dkd => 0,
        :dpe => 0,
        :salientpole => 1,
        :dpu => 0,
        :addmt => 0,
        :xmdm => 0,
    )
)

primary_parameters_gen2 = Dict(
    Symbol("machine__", k) => v for (k, v) in
    (
        :S_b => 100e6,
        :Sn => 192e6,
        :V_b => 18,
        :Vn => 18,
        :ω_b => 2π*60,
        :H => 3.921568,
        :D => 0,
        :R_s => 0.005,
        :X_rld => 0,
        :X_rlq => 0,
        :X_d => 1.719936,
        :X_q => 1.65984,
        :X′_d => 0.230016,
        :X′_q => 0.378048,
        :X″_d => 0.2,
        :X″_q => 0.2,
        :X_ls => 0.100032,
        :T′_d0 => 6,
        :T″_d0 => 0.0575,
        :T′_q0 => 0.535,
        :T″_q0 => 0.0945,
        :cosn => 0.85,
        :dkd => 0,
        :dpe => 0,
        :salientpole => 0,
        :dpu => 0,
        :addmt => 0,
        :xmdm => 0,
    )
)
avr_parameters_gen2 = Dict(
    Symbol("avr__", k) => v for (k, v) in
    (
        :Ka => 25,
        :Ke => -0.044,
        :Kf => 0.0805,
        :Ta => 0.2,
        :Tf => 0.35,
        :Te => 0.5,
        :Tr => 0.06,
        :vr_min => -1,
        :vr_max => 1,
        :anti_windup => false,
        :A => 0.0016,
        :B => 1.465,
    )
)

primary_parameters_gen3 = Dict(
    Symbol("machine__", k) => v for (k, v) in
    (
        :S_b => 100e6,
        :Sn => 128e6,
        :V_b => 13.8,
        :Vn => 13.8,
        :ω_b => 2π*60,
        :H => 2.766544,
        :D => 0,
        :R_s => 0.0001,
        :X_rld => 0,
        :X_rlq => 0,
        :X_d => 1.68,
        :X_q => 1.609984,
        :X′_d => 0.232064,
        :X′_q => 0.32,
        :X″_d => 0.2,
        :X″_q => 0.2,
        :X_ls => 0.094976,
        :T′_d0 => 5.89,
        :T″_d0 => 0.0575,
        :T′_q0 => 0.6,
        :T″_q0 => 0.08,
        :cosn => 0.85,
        :dkd => 0,
        :dpe => 0,
        :salientpole => 0,
        :dpu => 0,
        :addmt => 0,
        :xmdm => 0,
    )
)


@named mtkbus1 = StandardBus(; primary_parameters_gen1...)
@named mtkbus2 = StandardBus_AVR(; primary_parameters_gen2..., avr_parameters_gen2...)
@named mtkbus3 = StandardBus(; primary_parameters_gen3...)
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

#generate lines
t27_data = (;
    S_b = 100000000,
    S_n = 200000000,
    #V_OS = 230000,
    #V_US = 180000,
    X_pu_PF = 0.1250,
    R_pu_PF = 0
)
t27_params = trafoparams_pu(; t27_data...)

l45_data = (;
    S_b = 100000000,
    V_b = 230000,
    #ω_b = 60 * 2 * π,
    line_length = 1,
    R_l_perkm = 5.29,
    X_l_perkm = 44.965,
    B_l_perkm = 332.7 * 10^(-6)
)
l45_params = lineparams_pu(; l45_data...)

l57_data = (;
    S_b = 100000000,
    V_b = 230000,
    #ω_b = 60 * 2 * π,
    line_length = 1,
    R_l_perkm = 16.928,
    X_l_perkm = 85.169,
    B_l_perkm = 578.45 * 10^(-6),
    R_f_Ω=20, #not used
    X_f_Ω=30  #not used
)
l57_params = lineparams_pu(; l57_data...)

@named l46 = Line(piline(; R=0.0170, X=0.0920, B=0.1580), src=4, dst=6) #here in pu
@named l69 = Line(piline(; R=0.0390, X=0.1700, B=0.3580), src=6, dst=9)
@named l78 = Line(piline(; R=0.0085, X=0.0720, B=0.1490), src=7, dst=8)
@named l89 = Line(piline(; R=0.0119, X=0.1008, B=0.2090), src=8, dst=9)
@named t14 = Line(transformer(; R=0, X=0.0576), src=1, dst=4)
@named t27 = Line(transformer(; R=0, X=0.0625), src=2, dst=7)
@named t39 = Line(transformer(; R=0, X=0.0586), src=3, dst=9)
@named l57 = Line(piline_shortcircuit(; R=l57_params.R, X=l57_params.X, B=l57_params.B, pos=0.99, faultimp=0), src=5, dst=7) #faultimp=1, G_fault=l57_params.G_f, B_fault=l57_params.B_f
@named t27 = Line(transformer(; R=t27_params.R, X=t27_params.X), src=2, dst=7)
@named l45 = Line(piline(R=l45_params.R, X=l45_params.X, B=l45_params.B), src=4, dst=5)


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
    if integrator.t == 1
        @info "Short circuit on line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊shortcircuit] = 1
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_shortcircuit = PresetTimeCallback([1.0], affect1!)

affect2! = (integrator) -> begin
    if integrator.t == 1.05
        @info "Deactivate line 57 at t = $(integrator.t)"
        p = NWParameter(integrator)
        p.e[6, :pibranch₊active] = 0
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    else
        error("Should not be reached.")
    end
end
cb_deactivate = PresetTimeCallback([1.05], affect2!)

cb_set = CallbackSet(cb_shortcircuit, cb_deactivate)
prob = ODEProblem(nw, uflat(u0), (0,50), copy(pflat(u0)) ; callback=cb_set)
sol = solve(prob, Rodas5P(), dtmax=0.0001);
nothing

break # stop execution of script here

inspect(sol)

#Plot results
#### Voltage Magnitude
ref_bus57 = CSV.read("test/PFvalidation/PFdata/Bus5-7_standardModelPF_avrAmplidyne.csv", DataFrame; header=2, decimal=',')
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude (Power Factory Standard Model)")
ts = range(sol.t[begin],sol.t[end],length=10000)
umag5 = sqrt.(sol(ts; idxs=VIndex(5, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(5, :busbar₊u_i)).^2)
umag7 = sqrt.(sol(ts; idxs=VIndex(7, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(7, :busbar₊u_i)).^2)
lines!(ax, ts, umag5.u; label="Bus5")
lines!(ax, ref_bus57."Zeitpunkt in s", ref_bus57."u1, Betrag in p.u._1", color=Cycled(1), linestyle=:dash, label="Bus 5 ref")
lines!(ax, ts, umag7.u; label="Bus7")
lines!(ax, ref_bus57."Zeitpunkt in s", ref_bus57."u1, Betrag in p.u.", color=Cycled(2), linestyle=:dash, label="Bus 7 ref")
axislegend(ax; position=:rb)
xlims!(ax, 0, 5)
fig

### magnitude at generator bus
ref_bus = CSV.read("test/PFvalidation/PFdata/bus_voltmag_avrAmplidyne.csv", DataFrame; header=2, decimal=',') #_faultimp
fig = Figure();
ax = Axis(fig[1, 1]; title="Bus voltage magnitude")
ts = range(sol.t[begin],sol.t[end],length=10000)
umag1 = sqrt.(sol(ts; idxs=VIndex(1, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(1, :busbar₊u_i)).^2)
umag2 = sqrt.(sol(ts; idxs=VIndex(2, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(2, :busbar₊u_i)).^2)
umag3 = sqrt.(sol(ts; idxs=VIndex(3, :busbar₊u_r)).^2 + sol(ts; idxs=VIndex(3, :busbar₊u_i)).^2)
lines!(ax, ts, umag1.u; label="Bus 1")
lines!(ax, ref_bus."Zeitpunkt in s", ref_bus."Bus 1", color=Cycled(1), linestyle=:dash, label="Bus 1 ref")
lines!(ax, ts, umag2.u; label="Bus 2")
lines!(ax, ref_bus."Zeitpunkt in s", ref_bus."Bus 2", color=Cycled(2), linestyle=:dash, label="Bus 2 ref")
lines!(ax, ts, umag3.u; label="Bus 3")
lines!(ax, ref_bus."Zeitpunkt in s", ref_bus."Bus 3", color=Cycled(3), linestyle=:dash, label="Bus 3 ref")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 1.5)
fig


# Bus 2
ref_gen2 = CSV.read("test/PFvalidation/PFdata/gen2_data_avrAmplidyne.csv", DataFrame; header=2, decimal=',') #_faultimp
#### id and iq generator
fig = Figure();
ax = Axis(fig[1, 1]; title="stator current gen 2")
ts = range(sol.t[begin],sol.t[end],length=10000)
id = sol(ts; idxs=VIndex(2, :machine₊I_d))
iq = sol(ts; idxs=VIndex(2, :machine₊I_q))
lines!(ax, ts, id.u; label="i_d")
lines!(ax, ref_gen2."Zeitpunkt in s", ref_gen2."Ständerstrom, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="i_d ref")
lines!(ax, ts, iq.u; label="i_q")
lines!(ax, ref_gen2."Zeitpunkt in s", ref_gen2."Ständerstrom, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="i_q ref")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 10)
fig

#### ud and uq:
fig = Figure();
ax = Axis(fig[1, 1]; title="voltage at generator 2")
ts = range(sol.t[begin],sol.t[end],length=10000)
vd = sol(ts; idxs=VIndex(2, :machine₊V_d))
vq = sol(ts; idxs=VIndex(2, :machine₊V_q))
lines!(ax, ts, vd.u; label="u_d")
lines!(ax, ref_gen2."Zeitpunkt in s", ref_gen2."Spannung, d-Achse in p.u.", color=Cycled(1), linestyle=:dash, label="u_d ref")
lines!(ax, ts, vq.u; label="u_q")
lines!(ax, ref_gen2."Zeitpunkt in s", ref_gen2."Spannung, q-Achse in p.u.", color=Cycled(2), linestyle=:dash, label="u_q ref")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 2)
fig

#AVR data
ref_avr = CSV.read("test/PFvalidation/PFdata/Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';') #_faultimp
#vr in OpPoDyn
fig = Figure();
ax = Axis(fig[1, 1]; title="vr")
ts = range(sol.t[begin],sol.t[end],length=10000)
vr = sol(ts; idxs=VIndex(2, :avr₊vr))
vfout = sol(ts; idxs=VIndex(2, :avr₊vfout))
lines!(ax, ts, vr.u; label="vr")
lines!(ax, ref_avr."Zeitpunkt in s", ref_avr."o1", color=Cycled(1), linestyle=:dash, label="vr in PowerFactory")
#lines!(ax, ts, vfout.u; label="vfout")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 5)
fig

#v_mag.u
fig = Figure();
ax = Axis(fig[1, 1]; title="v_mag.u")
ts = range(sol.t[begin],sol.t[end],length=10000)
v_mag = sol(ts; idxs=VIndex(2, :machine₊v_mag))
lines!(ax, ts, v_mag.u; label="v_mag.u in OpPoDyn")
lines!(ax, ref_avr."Zeitpunkt in s", ref_avr."u", color=Cycled(1), linestyle=:dash, label="u in PowerFactory")
axislegend(ax; position=:rb)
xlims!(ax, 0.9, 5)
fig

#vref passt
ref_avr.summe = ref_avr."upss" .+ ref_avr."voel" .+ ref_avr."vuel" .+ ref_avr."avrref"
fig = Figure();
ax = Axis(fig[1, 1]; title="v_ref")
ts = range(sol.t[begin],sol.t[end],length=10000)
vref = sol(ts; idxs=VIndex(2, :avr₊vref)) ##möglich, da eig _vref, aber vref_input=false
lines!(ax, ts, vref.u; label="v_ref in OpPoDyn")
lines!(ax, ref_avr."Zeitpunkt in s", ref_avr.summe, color=Cycled(1), linestyle=:dash, label="v_ref aus PowerFactory")
axislegend(ax; position=:rt)
xlims!(ax, 0, 5)
fig

#output
ref_avr = CSV.read("test/PFvalidation/PFdata/Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';') #_faultimp
fig = Figure();
ax = Axis(fig[1, 1]; title="vfout")
ts = range(sol.t[begin],sol.t[end],length=10000)
vfout = sol(ts; idxs=VIndex(2, :avr₊vfout))
lines!(ax, ts, vfout.u; label="vfout")
lines!(ax, ref_avr."Zeitpunkt in s", ref_avr."uerrs", color=Cycled(1), linestyle=:dash, label="vfout in PowerFactory")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 10)
fig

#v_fb
fig = Figure();
ax = Axis(fig[1, 1]; title="v_fb")
ts = range(sol.t[begin],sol.t[end],length=10000)
vfout = sol(ts; idxs=VIndex(2, :avr₊v_fb))
lines!(ax, ts, vfout.u; label="v_fb")
lines!(ax, ref_avr."Zeitpunkt in s", ref_avr."vf", color=Cycled(1), linestyle=:dash, label="vf in PowerFactory")
axislegend(ax; position=:rt)
xlims!(ax, 0.9, 5)
fig

#test deviation from Power Factory (caution: only works with simulation time horizon (0,50))
function states_deviation(i, rmssym, ndsym)
    ref = CSV.read("test/PFvalidation/PFdata/Gen2_standardModelPF_avrdata.csv", DataFrame; header=2, decimal=',', delim=';') #_faultimp
    ref_t = ref[!, "Zeitpunkt in s"]  # Zeitwerte aus CSV
    ref_v = ref[!, rmssym] # Referenzwerte aus CSV
    sim_v = sol(ref_t, idxs=VIndex(i, ndsym)).u  # Simulation an den gleichen Zeitpunkten auswerten
    sum((sim_v - ref_v).^2) / length(ref_v)  # Mittelwert des quadratischen Fehlers
end

@testset "generator state deviation" begin
    @test states_deviation(2, :vf, :avr₊v_fb) < 1e-5
    @test states_deviation(2, :o1, :avr₊vr) < 1e-3
    @test states_deviation(2, :uerrs, :avr₊vfout) < 1e-4
    @test states_deviation(2, :u, :machine₊v_mag) < 1e-5
end