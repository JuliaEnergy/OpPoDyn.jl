module OpPoDyn  # v2

using Reexport: Reexport, @reexport
@reexport using NetworkDynamics
using NetworkDynamics: SymbolicView

@reexport using PowerDynamics

using SciMLBase: SciMLBase, solve
using ForwardDiff: ForwardDiff
using ArgCheck: @argcheck
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII

using ModelingToolkit: ModelingToolkit, @connector, @mtkmodel, @variables, @named,
                       System, connect, getname, unknowns, get_name, get_iv, get_systems,
                       get_gui_metadata, t_nounits as t, Equation,
                       defaults, parameters, iscomplete, rename, simplify, unwrap,
                       get_eqs, get_observed, get_defaults, get_schedule,
                       get_connector_type, get_preface, get_initializesystem,
                       get_continuous_events, get_discrete_events, get_parameter_dependencies,
                       get_tspan, get_guesses,
                       structural_simplify
using ModelingToolkit: @unpack, Num, System # needed for @mtkmodel?
using Symbolics: Symbolics, Symbolic, iscall, fixpoint_sub
using SciMLBase: SciMLBase, solve

include("Library/Library.jl")

export pin_parameters
include("pin_parameters.jl")

end
