"""
Circuit.jl

Series circuit of TRIOMA components (Component, BreedingBlanket, GLC).
Supports open and closed-loop (iterative) solving.
"""

module CircuitModule

using ..TriomaTypes
using ..PipeSubclasses
using ..PAVModule
import ..PAVModule: get_inventory!, estimate_cost!
using ..BreedingBlanketModule
import ..BreedingBlanketModule: connect_to_component!
using ..GasLiquidContactor

export Circuit,
       add_component!,
       solve_circuit!,
       get_eff_circuit!,
       get_gains_and_losses!,
       get_circuit_pumping_power!,
       inspect_circuit

# ---------------------------------------------------------------------------
# Union type alias for the three valid component types
# ---------------------------------------------------------------------------
const AnyComponent = Union{Component, BreedingBlanket, GLC}

# ---------------------------------------------------------------------------
# Circuit struct
# ---------------------------------------------------------------------------

"""
Series circuit of TRIOMA components.

Fields:
- `components`    : ordered `Vector` of `Component`, `BreedingBlanket`, or `GLC` objects
- `closed`        : if `true`, the circuit forms a closed loop and `solve_circuit!`
                    iterates until the inlet concentration of the first component
                    converges with the outlet of the last
- `eff`           : overall extraction efficiency (set by `get_eff_circuit!` or `get_gains_and_losses!`)
- `extraction_perc` : fraction of tritium extracted (set by `get_gains_and_losses!`)
- `loss_perc`     : fraction of tritium lost (set by `get_gains_and_losses!`)
- `inv`           : total tritium inventory [mol] (set by `get_inventory!`)
- `pumping_power` : total pumping power [W] (set by `get_circuit_pumping_power!`)
- `cost`          : total estimated cost (set by `estimate_cost!`)
"""
mutable struct Circuit <: TriomaClass
    components    ::Vector{AnyComponent}
    closed        ::Bool
    eff           ::Union{Float64, Nothing}
    extraction_perc::Union{Float64, Nothing}
    loss_perc     ::Union{Float64, Nothing}
    inv           ::Union{Float64, Nothing}
    pumping_power ::Union{Float64, Nothing}
    cost          ::Union{Float64, Nothing}
end

"""
    Circuit(; components=nothing, closed=false)

Construct a `Circuit`. `components` may be a `Vector` of `Component`,
`BreedingBlanket`, `GLC`, or nested `Circuit` objects (nested circuits are
flattened into a single list). Passing an invalid type throws `ArgumentError`.
"""
function Circuit(;
    components::Union{Vector, Nothing}=nothing,
    closed::Bool=false)

    vec = AnyComponent[]
    if components !== nothing
        for element in components
            if isa(element, AnyComponent)
                push!(vec, element)
            elseif isa(element, Circuit)
                append!(vec, element.components)
            else
                throw(ArgumentError(
                    "Invalid component type: $(typeof(element)). " *
                    "Expected Component, BreedingBlanket, GLC, or Circuit."))
            end
        end
    end
    return Circuit(vec, closed, nothing, nothing, nothing, nothing, nothing, nothing)
end

# ---------------------------------------------------------------------------
# connect_to_component! — the key wiring primitive
# ---------------------------------------------------------------------------

"""
    BreedingBlanketModule.connect_to_component!(src, dst)

Set `dst.c_in` equal to `src.c_out`. Works for any TRIOMA component pair.
This is the generic fallback; typed methods for `BreedingBlanket` and `GLC`
are defined in their own modules.
"""
function connect_to_component!(src, dst)
    dst === nothing && throw(ArgumentError("destination component cannot be nothing"))
    update_attribute!(dst, "c_in", src.c_out)
end

# ---------------------------------------------------------------------------
# add_component!
# ---------------------------------------------------------------------------

"""
    add_component!(circuit::Circuit, component)

Append a component to the circuit. If `component` is itself a `Circuit`,
its contents are flattened in.
"""
function add_component!(circuit::Circuit, component)
    if isa(component, Circuit)
        append!(circuit.components, component.components)
    elseif isa(component, AnyComponent)
        push!(circuit.components, component)
    else
        throw(ArgumentError(
            "Invalid component type: $(typeof(component))."))
    end
end

# ---------------------------------------------------------------------------
# solve_circuit!
# ---------------------------------------------------------------------------

"""
    solve_circuit!(circuit::Circuit; tol=1e-6)

March through the circuit in order, computing each component's outlet
concentration and propagating it to the next component's inlet.

- `Component`       : analytical efficiency → `outlet_c_comp!`
- `BreedingBlanket` : `get_cout!`
- `GLC`             : `get_c_out!`

If `circuit.closed == true`, the loop repeats until the relative change in
the first component's inlet concentration is below `tol`.

Note: the Python source warns when more than one `BreedingBlanket` is present;
this translation preserves that warning.
"""
function solve_circuit!(circuit::Circuit; tol::Float64=1e-6)
    comps  = circuit.components
    n      = length(comps)

    bb_count = count(c -> isa(c, BreedingBlanket), comps)
    bb_count > 1 && println("Warning: there are more than one BreedingBlanket in the circuit!")

    converged = false
    while !converged
        for (i, comp) in enumerate(comps)
            if isa(comp, GLC)
                get_c_out!(comp)
                i < n && connect_to_component!(comp, comps[i + 1])

            elseif isa(comp, Component)
                use_analytical_efficiency!(comp; p_out=comp.p_out)
                outlet_c_comp!(comp)
                i < n && connect_to_component!(comp, comps[i + 1])

            elseif isa(comp, BreedingBlanket)
                get_cout!(comp)
                i < n && connect_to_component!(comp, comps[i + 1])
            end
        end

        if circuit.closed
            err = abs(comps[1].c_in - comps[end].c_out) / comps[1].c_in
            converged = err < tol
        else
            converged = true
        end

        # Feed last outlet back to first inlet for next closed-loop iteration
        connect_to_component!(comps[end], comps[1])
    end
end

# ---------------------------------------------------------------------------
# get_eff_circuit!
# ---------------------------------------------------------------------------

"""
    get_eff_circuit!(circuit::Circuit)

Run one open-loop pass through the circuit (after rotating so the
`BreedingBlanket` is first) and compute the overall extraction efficiency
relative to the concentration entering the first non-BB component.

Stores the result in `circuit.eff`.
"""
function get_eff_circuit!(circuit::Circuit)
    comps = circuit.components
    n     = length(comps)

    # Rotate so BreedingBlanket is first
    bb_idx = findfirst(c -> isa(c, BreedingBlanket), comps)
    if bb_idx !== nothing
        circuit.components = vcat(comps[bb_idx:bb_idx],
                                  comps[1:bb_idx-1],
                                  comps[bb_idx+1:end])
        comps = circuit.components
    end

    for (i, comp) in enumerate(comps)
        if isa(comp, GLC)
            get_c_out!(comp)
            i < n && connect_to_component!(comp, comps[i + 1])
        end
        if isa(comp, Component)
            use_analytical_efficiency!(comp; p_out=comp.p_out)
            outlet_c_comp!(comp)
        end
        i < n && connect_to_component!(comp, comps[i + 1])
    end

    # Efficiency is measured from the first post-BB component
    circuit.eff = (comps[2].c_in - comps[end].c_out) / comps[2].c_in
end

# ---------------------------------------------------------------------------
# get_gains_and_losses!
# ---------------------------------------------------------------------------

"""
    get_gains_and_losses!(circuit::Circuit)

Tally the concentration drop across each `Component`, splitting extractors
(`loss == false`) from loss pipes (`loss == true`). Normalises by the
concentration entering the first post-BB component and by the overall
circuit efficiency.

Stores `circuit.extraction_perc`, `circuit.loss_perc`, and `circuit.eff`.
"""
function get_gains_and_losses!(circuit::Circuit)
    comps  = circuit.components
    n      = length(comps)

    gains  = sum(c.c_in - c.c_out for c in comps if isa(c, Component) && !c.loss; init=0.0)
    losses = sum(c.c_in - c.c_out for c in comps if isa(c, Component) &&  c.loss; init=0.0)

    bb_count = count(c -> isa(c, BreedingBlanket), comps)
    bb_count > 1 && println("Warning: there are more than one BreedingBlanket!")
    bb_idx = findfirst(c -> isa(c, BreedingBlanket), comps)

    # Determine reference inlet/outlet, skipping the BreedingBlanket
    if bb_idx === nothing
        ref_in, ref_out = comps[1].c_in,       comps[end].c_out
    elseif bb_idx == 1
        ref_in, ref_out = comps[2].c_in,       comps[end].c_out
    elseif bb_idx == n
        ref_in, ref_out = comps[1].c_in,       comps[n-1].c_out
    else
        ref_in, ref_out = comps[bb_idx+1].c_in, comps[bb_idx-1].c_out
    end

    eff_circuit             = (ref_in - ref_out) / ref_in
    circuit.extraction_perc = gains  / ref_in / eff_circuit
    circuit.loss_perc       = losses / ref_in / eff_circuit
    circuit.eff             = eff_circuit
end

# ---------------------------------------------------------------------------
# get_inventory!
# ---------------------------------------------------------------------------

"""
    PAVModule.get_inventory!(circuit::Circuit; flag_an=true)

Sum the tritium inventory across all `Component` objects in the circuit.
Stores the total in `circuit.inv` [mol].
"""
function get_inventory!(circuit::Circuit; flag_an::Bool=true)
    total = 0.0
    for comp in circuit.components
        if isa(comp, Component)
            get_inventory!(comp; flag_an=flag_an)
            total += comp.inv
        end
    end
    circuit.inv = total
end

# ---------------------------------------------------------------------------
# get_circuit_pumping_power!
# ---------------------------------------------------------------------------

"""
    get_circuit_pumping_power!(circuit::Circuit) -> Float64

Sum the pumping power of all `Component` objects in the circuit.
Stores the result in `circuit.pumping_power` [W] and returns it.
"""
function get_circuit_pumping_power!(circuit::Circuit)::Float64
    total = 0.0
    for comp in circuit.components
        if isa(comp, Component)
            get_pressure_drop!(comp)
            get_pumping_power!(comp)
            total += comp.pumping_power
        end
    end
    circuit.pumping_power = total
    return circuit.pumping_power
end

# ---------------------------------------------------------------------------
# inspect_circuit
# ---------------------------------------------------------------------------

"""
    inspect_circuit(circuit::Circuit; name=nothing)

Print all fields of each component. If `name` is provided, only print the
component whose `name` field matches.
"""
function inspect_circuit(circuit::Circuit; name::Union{String, Nothing}=nothing)
    for comp in circuit.components
        if name === nothing || (hasproperty(comp, :name) && getfield(comp, :name) == name)
            inspect(comp)
        end
    end
end

# ---------------------------------------------------------------------------
# estimate_cost!
# ---------------------------------------------------------------------------

"""
    PAVModule.estimate_cost!(circuit::Circuit; metal_costs, fluid_costs) -> Float64

Estimate total fabrication cost across all `Component` objects in the circuit.
`metal_costs` and `fluid_costs` must be indexable vectors aligned with
`circuit.components` (non-`Component` entries are skipped).

Stores total in `circuit.cost` and returns it.
"""
function estimate_cost!(circuit::Circuit;
                        metal_costs::Vector{Float64},
                        fluid_costs::Vector{Float64})::Float64
    total = 0.0
    for (i, comp) in enumerate(circuit.components)
        if isa(comp, Component)
            estimate_cost!(comp; metal_cost=metal_costs[i], fluid_cost=fluid_costs[i])
            total += comp.cost
        end
    end
    circuit.cost = total
    return circuit.cost
end

end # module
