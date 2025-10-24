module Library

using ArgCheck: @argcheck
using ModelingToolkit: ModelingToolkit, @named, @mtkmodel, @variables, @parameters, simplify, connect,
                       t_nounits as t, D_nounits as Dt
using ModelingToolkit: @unpack, Equation, Num, System # needed for @mtkmodel?
using ModelingToolkitStandardLibrary.Blocks #: RealInput, RealOutput
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve

using PowerDynamics
using PowerDynamics.Library

@mtkmodel SystemBase begin
    @parameters begin
        SnRef = 100, [description="System base"]
        fNom = 50, [description="AC system frequency"]
        ωNom = 2 * π * fNom, [description="System angular frequency"]
        ωRef0Pu = 1, [description="Reference for system angular frequency (pu base ωNom)"]
        ω0Pu = 1, [description="System angular frequency (pu base ωNom)"]
   end
end

####
#### Machine Models
####
export DynawoMachine
include("Machines/DynawoMachine.jl")

export IPSLPSATOrder4
include("Machines/IPSLPSAT.jl")

export ClassicalMachine_powerfactory
include("Machines/ClassicalMachine_powerfactory.jl")

export ClassicalMachine_pf
include("Machines/ClassicalMachine_pf.jl")

export StandardModel_pf_testneu
include("Machines/StandardModel_pf_testneu.jl")

export StandardModel_pf
include("Machines/StandardModel_pf.jl")

####
#### Line Models
####
export DynawoPiLine
include("Branches/DynawoPiLine.jl")

export DynawoFixedRatioTransformer
include("Transformers/DynawoFixedRatioTransformer.jl")

####
#### Load Models
####
#export PQLoad, VoltageDependentLoad, ConstantYLoad, ZIPLoad
#include("Loads/PQLoad.jl")

####
#### WECC modules
####
export limiter, lowlimit, uplimit, deadband, LVPLogic, VDL
include("WECC-models/functions.jl")

export regc_a
include("WECC-models/regc.jl")

export reec_a, reec_b, reec_c
include("WECC-models/reec.jl")

export repc_a
include("WECC-models/repc.jl")

export WECC_large_PV, WECC_BESS, WECC_WT_4B, WTDTA1
include("WECC-models/plantmodels.jl")

export OpenIPSL_RePSSE, ref_rms_error
include("OpenIPSL/test/OpenIPSL_testenvRenewablePSSE.jl")

end

