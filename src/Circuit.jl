"""
    CircuitModule

Series circuit of fuel-cycle components with solve, efficiency, and inventory methods.
"""
module CircuitModule

using ..TriomaModule: TriomaClass, update_attribute!
using ..PAVComponent: Component, use_analytical_efficiency!, outlet_c_comp!, get_inventory!
using ..BreedingBlanketModule: BreedingBlanket, get_cout!
using ..GasLiquidContactor: GLC, get_c_out!

export Circuit, add_component!, get_eff_circuit!, get_gains_and_losses!,
       solve_circuit!, get_circuit_inventory!, get_circuit_pumping_power!,
       inspect_circuit, estimate_circuit_cost

"""Union of all component types accepted by Circuit."""
const AnyComponent = Union{Component, BreedingBlanket, GLC}

# Helper: forward connect_to_component (sets downstream c_in = upstream c_out)
function _connect!(upstream::AnyComponent, downstream::AnyComponent)
    update_attribute!(downstream, :c_in, upstream.c_out)
end

"""
    Circuit

A series chain of fuel-cycle components.

Fields
- `components` : ordered `Vector` of `Component`, `BreedingBlanket`, or `GLC` objects
- `closed`     : if `true`, the circuit is a closed loop (solved iteratively)
- `eff`        : overall tritium extraction efficiency (set by `get_eff_circuit!`)
- `inv`        : total tritium inventory [mol] (set by `get_circuit_inventory!`)
- `pumping_power` : total pumping power [W]
- `cost`       : total cost estimate [USD]
- `extraction_perc`, `loss_perc` : fractions from `get_gains_and_losses!`
"""
mutable struct Circuit <: TriomaClass
    components::Vector{AnyComponent}
    closed::Bool
    eff::Union{Float64,Nothing}
    inv::Union{Float64,Nothing}
    pumping_power::Union{Float64,Nothing}
    cost::Union{Float64,Nothing}
    extraction_perc::Union{Float64,Nothing}
    loss_perc::Union{Float64,Nothing}
end

function Circuit(; components=nothing, closed::Bool=false)
    vec_comp = AnyComponent[]
    if components !== nothing
        for el in components
            if el isa AnyComponent
                push!(vec_comp, el)
            elseif el isa Circuit
                append!(vec_comp, el.components)
            else
                throw(ArgumentError("Invalid component type: $(typeof(el))"))
            end
        end
    end
    Circuit(vec_comp, closed, nothing, nothing, nothing, nothing, nothing, nothing)
end

"""Append one or more components (or flatten a nested Circuit) to the circuit."""
function add_component!(circ::Circuit, component::Union{AnyComponent,Circuit})
    if component isa Circuit
        append!(circ.components, component.components)
    else
        push!(circ.components, component)
    end
end

"""
    get_eff_circuit!(circ)

Propagate concentrations forward through all components and compute overall
extraction efficiency.  A `BreedingBlanket`, if present, is moved to position 1.
"""
function get_eff_circuit!(circ::Circuit)
    # Move BreedingBlanket to the front if present
    bb_idx = findfirst(c -> c isa BreedingBlanket, circ.components)
    if bb_idx !== nothing
        bb = circ.components[bb_idx]
        deleteat!(circ.components, bb_idx)
        pushfirst!(circ.components, bb)
    end

    for (i, comp) in enumerate(circ.components)
        if comp isa GLC
            get_c_out!(comp)
            i < length(circ.components) && _connect!(comp, circ.components[i+1])
        elseif comp isa Component
            use_analytical_efficiency!(comp; p_out=comp.p_out)
            outlet_c_comp!(comp)
            i < length(circ.components) && _connect!(comp, circ.components[i+1])
        end
    end

    circ.eff = (circ.components[2].c_in - circ.components[end].c_out) /
               circ.components[2].c_in
end

"""
    get_gains_and_losses!(circ)

Decompose circuit efficiency into extraction (gains) and loss fractions.
"""
function get_gains_and_losses!(circ::Circuit)
    gains = 0.0; losses = 0.0
    bb_idx = nothing

    for (i, comp) in enumerate(circ.components)
        comp isa BreedingBlanket && (bb_idx = i)
        if comp isa Component
            diff = comp.c_in - comp.c_out
            comp.loss ? (losses += diff) : (gains += diff)
        end
    end

    bb_idx === nothing && error("No BreedingBlanket found in circuit")

    if bb_idx == 1
        c_ref = circ.components[2].c_in
        eff   = (circ.components[2].c_in - circ.components[end].c_out) / c_ref
    elseif bb_idx == length(circ.components)
        c_ref = circ.components[1].c_in
        eff   = (circ.components[1].c_in - circ.components[bb_idx-1].c_out) / c_ref
    else
        c_ref = circ.components[bb_idx+1].c_in
        eff   = (circ.components[bb_idx+1].c_in - circ.components[bb_idx-1].c_out) / c_ref
    end

    circ.eff            = eff
    circ.extraction_perc = gains   / c_ref / eff
    circ.loss_perc       = losses  / c_ref / eff
end

"""
    solve_circuit!(circ; tol=1e-6)

Solve the circuit by propagating concentrations forward.  For closed loops
(`circ.closed = true`) iterates until the inlet concentration converges.
"""
function solve_circuit!(circ::Circuit; tol::Float64=1e-6)
    bb_idx = findfirst(c -> c isa BreedingBlanket, circ.components)

    flag = false
    while !flag
        for (i, comp) in enumerate(circ.components)
            if comp isa GLC
                get_c_out!(comp)
            elseif comp isa Component
                use_analytical_efficiency!(comp; p_out=comp.p_out)
                outlet_c_comp!(comp)
            elseif comp isa BreedingBlanket
                get_cout!(comp)
            end
            i < length(circ.components) && _connect!(comp, circ.components[i+1])
        end

        if circ.closed
            c0  = circ.components[1].c_in
            c_f = circ.components[end].c_out
            err = abs(c0 - c_f) / c0
            err < tol && (flag = true)
            _connect!(circ.components[end], circ.components[1])
        else
            flag = true
        end
    end
end

"""
    get_circuit_inventory!(circ; flag_an=true)

Compute the total tritium inventory [mol] across all `Component` members.
"""
function get_circuit_inventory!(circ::Circuit; flag_an::Bool=true)
    inv = 0.0
    for comp in circ.components
        comp isa Component && (get_inventory!(comp; flag_an=flag_an); inv += comp.inv)
    end
    circ.inv = inv
end

"""
    get_circuit_pumping_power!(circ) -> pumping_power [W]

Compute total pumping power across all `Component` members.
"""
function get_circuit_pumping_power!(circ::Circuit)
    using ..PAVComponent: get_pressure_drop!, get_pumping_power!
    power = 0.0
    for comp in circ.components
        if comp isa Component
            get_pressure_drop!(comp)
            get_pumping_power!(comp)
            power += comp.pumping_power
        end
    end
    circ.pumping_power = power
    return power
end

"""Print all (or named) components in the circuit."""
function inspect_circuit(circ::Circuit; name::Union{String,Nothing}=nothing)
    using ..TriomaModule: inspect
    for comp in circ.components
        if name === nothing || (hasproperty(comp, :name) && comp.name == name)
            inspect(comp)
        end
    end
end

"""
    estimate_circuit_cost(circ; metal_costs, fluid_costs) -> total_cost [USD]

`metal_costs` and `fluid_costs` are vectors ordered the same as `circ.components`.
"""
function estimate_circuit_cost(circ::Circuit; metal_costs::Vector, fluid_costs::Vector)
    using ..PAVComponent: estimate_cost
    total = 0.0
    for (i, comp) in enumerate(circ.components)
        if comp isa Component
            estimate_cost(comp; metal_cost=metal_costs[i], fluid_cost=fluid_costs[i])
            total += comp.cost
        end
    end
    circ.cost = total
    return total
end

end # module CircuitModule
