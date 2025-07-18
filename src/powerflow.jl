"""
    pfSlack(; V=missing, δ=missing, u_r=missing, u_i=missing)

Create a slack bus for power flow analysis.

A slack bus maintains constant voltage magnitude and phase angle (or real and imaginary voltage components).
Either provide voltage magnitude `V` and phase angle `δ`, or provide real and imaginary voltage components `u_r` and `u_i`.
"""
function pfSlack(; V=missing, δ=missing, u_r=missing, u_i=missing)
    if !ismissing(V) && ismissing(u_r) && ismissing(u_i)
        δ = ismissing(δ) ? 0 : δ
        @named slack = Library.VδConstraint(; V, δ)
        u_r = V * cos(δ)
        u_i = V * sin(δ)
    elseif ismissing(V) && ismissing(δ) && !ismissing(u_r) && !ismissing(u_i)
        @named slack = Library.UrUiConstraint(; u_r, u_i)
    else
        throw(ArgumentError("Either Provide V or δ, or u_r and u_i. But not both!"))
    end

    mtkbus = MTKBus(slack, name=:slackbus)
    b = Bus(mtkbus)
    set_voltage!(b, u_r + im * u_i)
    b
end

"""
    pfPV(; P, V)

Create a PV bus for power flow analysis.

A PV bus maintains constant active power injection and voltage magnitude.
The reactive power and voltage phase angle are determined by the power flow solution.
"""
function pfPV(; P, V)
    @named pv = Library.PVConstraint(; P, V)
    mtkbus = MTKBus(pv, name=:pvbus)
    b = Bus(mtkbus)
    set_voltage!(b; mag=V, arg=0)
    b
end

"""
    pfPQ(; P=0, Q=0)

Create a PQ bus for power flow analysis.

A PQ bus has specified active and reactive power injections.
The voltage magnitude and phase angle are determined by the power flow solution.
"""
function pfPQ(; P=0, Q=0)
    @named pq = Library.PQConstraint(; P, Q)
    mtkbus = MTKBus(pq, name=:pqbus)
    b = Bus(mtkbus)
    set_voltage!(b; mag=1, arg=0)
    b
end

"""
    ispfmodel(cf::NetworkDynamics.ComponentModel)

Check if a component model is suitable for power flow analysis.

A component model is considered a valid power flow model if it has no dynamics,
i.e., either no states or a zero mass matrix.

## Returns
- `true` if the component is suitable for power flow analysis
- `false` otherwise
"""
function ispfmodel(cf::NetworkDynamics.ComponentModel)
    # no states or mass matrix 0
    if NetworkDynamics.dim(cf) == 0 || cf.mass_matrix==LinearAlgebra.UniformScaling(0)
        return true
    end
    return false
end

"""
    powerflow_model(cf::NetworkDynamics.ComponentModel)

Extract or create a power flow component model from a dynamic component model.

1. If the component has `:pfmodel` metadata, use that model (after validation)
2. If the component is already a valid power flow model (i.e. no ODE, just constraints), return it as-is

## Returns
- A component model suitable for power flow analysis (no dynamics)

## Validation
The returned model must satisfy [`ispfmodel`](@ref) criteria:
- Either no states or zero mass matrix (no dynamics)

See also: [`ispfmodel`](@ref), [`pfSlack`](@ref), [`pfPV`](@ref), [`pfPQ`](@ref)
"""
function powerflow_model(cf::NetworkDynamics.ComponentModel)
    if has_metadata(cf, :pfmodel)
        pfm = get_metadata(cf, :pfmodel)
        if !ispfmodel(pfm)
            error("Provided :pfmodel for :$(cf.name) is no valid powerflow model!")
        end
        return pfm
    elseif ispfmodel(cf)
        return cf
    elseif cf isa VertexModel && has_default(cf, :busbar₊u_r) && has_default(cf, :busbar₊u_i)
        @warn "No powerflow model given for :$(cf.name), using slack with default voltage!"
        return pfSlack(u_r=get_default(cf, :busbar₊u_r), u_i=get_default(cf, :busbar₊u_i))
    end
    error("Cannot create PF component model from :$(cf.name)! Please proved :pfmodel metadata!")
end

"""
    powerflow_model(nw::Network)

Create a power flow network model from a dynamic network model.

This method applies [`powerflow_model`](@ref) to all vertex and edge components
in the network, creating a new network suitable for steady-state power flow analysis.

Returns a new `Network` with the same graph structure but power flow component models

See also: [`solve_powerflow`](@ref)
"""
function powerflow_model(nw::Network)
    g = nw.im.g
    vfs = powerflow_model.(nw.im.vertexm);
    efs = powerflow_model.(nw.im.edgem);
    Network(g, vfs, efs)
end

"""
    solve_powerflow(nw::Network;
                    pfnw = powerflow_model(nw),
                    pfs0 = NWState(nw),
                    verbose=true)

Solve the power flow equations for a given network.

Uses [`find_fixpoint`](@extref NetworkDynamics.find_fixpoint) from NetworkDynamics to solve the algebraic power flow equations.

## Parameters
- `nw`: The dynamic network model
- `pfnw`: The power flow network model (default: created from `nw`)
- `pfs0`: Initial state for the power flow calculation
- `verbose`: Whether to print the power flow solution

## Returns
- A `NWState` containing the solved power flow solution

See also [`initialize_from_pf`](@ref).
"""
function solve_powerflow(
    nw::Network;
    pfnw = powerflow_model(nw),
    pfs0 = NWState(nw),
    verbose=true
)
    pfnw.mass_matrix == LinearAlgebra.UniformScaling(0) || error("Powerflow model must have a mass matrix of 0!")

    uf = uflat(pfs0)
    pf = pflat(pfs0)
    pfs = find_fixpoint(pfnw, pfs0)
    verbose && show_powerflow(pfs)

    return pfs
end

initialize_from_pf_docstring = raw"""
    initialize_from_pf[!](
        nw::Network;
        verbose = true,
        subverbose = false,
        pfnw = powerflow_model(nw),
        pfs0 = NWState(pfnw),
        pfs = solve_powerflow(pfnw; pfs0, verbose),
        kwargs...
    )

Initialize a dynamic network model from a power flow solution.

This function performs a two-step initialization process:
1. Solve the power flow equations for the network
2. Use the power flow solution to initialize the dynamic model

There are two versions of this function: a mutating one (!-at the end of name) and a non-mutating version.
The mutating version uses `initialize_componentwise!` internally, the non-mutating one `initialize_componentwise`.
When the mutating version is used, `NWState(nw)` after initialization will return the same initialized
state again, as it is stored in the metadata.

## Parameters
- `nw`: The dynamic network model to initialize
- `verbose`: Whether to print information about the power flow solution (default: true)
- `subverbose`: Whether to print detailed information during component initialization (default: false). Can be Vector [VIndex(1), EIndex(3), ...] for selective output
- `pfnw`: Power flow network model (default: created from `nw` using `powerflow_model`)
- `pfs0`: Initial state for power flow calculation (default: created from `pfnw`)
- `pfs`: Power flow solution (default: calculated using `solve_powerflow`)
- Additional keyword arguments are passed to `initialize_componentwise[!]`

## Returns
- A fully initialized network state

See also: [`solve_powerflow`](@ref), [`initialize_componentwise`](@extref NetworkDynamics.initialize_componentwise), [`interface_values`](@extref NetworkDynamics.interface_values)
"""
@doc initialize_from_pf_docstring
initialize_from_pf(nw; kw...) = _init_from_pf(initialize_componentwise, nw; kw...)
@doc initialize_from_pf_docstring
initialize_from_pf!(nw; kw...) = _init_from_pf(initialize_componentwise!, nw; kw...)
function _init_from_pf(
    initf, nw;
    verbose = true,
    subverbose = false,
    pfnw = powerflow_model(nw),
    pfs0 = NWState(pfnw),
    pfs = solve_powerflow(pfnw; pfs0, verbose),
    kwargs...
)
    interface_vals = interface_values(pfs)
    pfinitconstraints = specialize_pfinitconstraints(nw, pfs)
    pfinitformulas = specialize_pfinitformulas(nw, pfs)
    initf(
        nw;
        default_overrides=interface_vals,
        additional_initconstraint = pfinitconstraints,
        additional_initformula = pfinitformulas,
        verbose, subverbose, kwargs...
    )
end

"""
    show_powerflow(s::NWState/Network)

Display power flow results in a tabular format.

Extract and format power flow solution data from a network state, showing bus-level information
including voltage magnitudes, phase angles, active power, and reactive power.
"""
function show_powerflow(s::NWState)
    NV = nv(extract_nw(s))
    dict = OrderedDict()
    dict["N"] = 1:NV
    dict["Bus Names"] = [cf.name for cf in extract_nw(s).im.vertexm]
    dict["vm [pu]"] = s[vidxs(1:NV, :busbar₊u_mag)]
    dict["varg [deg]"] = rad2deg.(s[vidxs(1:NV, :busbar₊u_arg)])
    dict["P [pu]"] = s[vidxs(1:NV, :busbar₊P)]
    dict["Q [pu]"] = s[vidxs(1:NV, :busbar₊Q)]

    DataFrame(dict)
end
show_powerflow(nw::Network) = show_powerflow(NWState(nw))
