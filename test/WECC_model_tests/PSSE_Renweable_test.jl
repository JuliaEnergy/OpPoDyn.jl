using PowerDynamics
PowerDynamics.load_pdtesting()
using Main.PowerDynamicsTesting
using OpPoDyn

using PowerDynamics.Library
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

using CSV
using DataFrames
using CairoMakie
#=
ref = CSV.read(
    joinpath(pkgdir(PowerDynamics),"test","OpenIPSL_test","GENROE","modelica_results.csv.gz"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
) =#

# bus 1 is provided from outside
PV_BUS = let
    ω_b = 2π*60

    # Powerflow results
    v_0 = 1.0
    angle_0 = 1.475 * π /180

    @named PV = OpPoDyn.Library.WECC_large_PV()
    busmodel = MTKBus(PV; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_RePSSE(PV_BUS);

BESS_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.0
    angle_0 = 1.475 * π /180 #in rad

    @named BESS = OpPoDyn.Library.WECC_BESS()
    busmodel = MTKBus(BESS; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_RePSSE(BESS_BUS);

WT4B_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.0
    angle_0 = 1.475 * π /180

    @named WT = OpPoDyn.Library.WECC_WT_4B()
    busmodel = MTKBus(WT; name=:GEN1)
    compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
end

sol = OpenIPSL_RePSSE(WT4B_BUS);
#=
## perform tests for all variables of interest
# Core machine variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊w), "gENROE.w") < 1e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊delta), "gENROE.delta") < 1e-3
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊P), "gENROE.P") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊Q), "gENROE.Q") < 5e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊Vt), "gENROE.Vt") < 2e-5

# State variables
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊Epd), "gENROE.Epd") < 3e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊Epq), "gENROE.Epq") < 6e-5
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊PSIkd), "gENROE.PSIkd") < 2e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊PSIkq), "gENROE.PSIkq") < 4e-4

# Field and torque
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊XadIfd), "gENROE.XadIfd") < 6e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊Te), "gENROE.Te") < 5e-4

# Current and voltage components
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊id), "gENROE.id") < 5e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊iq), "gENROE.iq") < 3e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊ud), "gENROE.ud") < 4e-4
@test ref_rms_error(sol, ref, VIndex(:GEN1, :genroe₊uq), "gENROE.uq") < 3e-4
=#