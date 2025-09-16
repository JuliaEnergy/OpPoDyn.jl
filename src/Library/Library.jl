module Library

using ArgCheck: @argcheck
using ..OpPoDyn: Terminal, BusBase, Ibase
using ModelingToolkit: ModelingToolkit, @named, @mtkmodel, @variables, @parameters, simplify,
                       t_nounits as t, D_nounits as Dt
using ModelingToolkit: @unpack, Equation, Num, System # needed for @mtkmodel?
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using NonlinearSolve: NonlinearProblem
using SciMLBase: SciMLBase, solve

using PowerDynamics.Library: BusBase

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

end
