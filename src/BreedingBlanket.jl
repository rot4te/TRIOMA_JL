"""
    BreedingBlanketModule

Breeding blanket component for tritium generation and outlet-concentration calculation.
"""
module BreedingBlanketModule

using ..TriomaModule: TriomaClass, update_attribute!
using ..PipeSubclasses: Fluid

export BreedingBlanket, get_flowrate!, get_cout!

const N_A              = 6.02214076e23   # Avogadro's number [1/mol]
const EV_TO_J          = 1.602176634e-19 # eV → J
const REACTION_ENERGY  = 17.6e6         # D-T fusion energy [eV]

"""
    BreedingBlanket

Represents a breeding blanket in the outer fuel cycle.

Fields
- `c_in`      : inlet tritium concentration [mol/m³]
- `Q`         : thermal power [W]
- `TBR`       : tritium breeding ratio [-]
- `T_out`     : fluid outlet temperature [K]
- `T_in`      : fluid inlet temperature [K]
- `fluid`     : `Fluid` object (for cp and rho)
- `name`      : optional label
- `m_coolant` : coolant mass flow rate [kg/s] (computed by `get_flowrate!`)
- `c_out`     : outlet tritium concentration [mol/m³] (set by `get_cout!`)
"""
mutable struct BreedingBlanket <: TriomaClass
    c_in::Union{Float64,Nothing}
    Q::Union{Float64,Nothing}
    TBR::Union{Float64,Nothing}
    T_out::Union{Float64,Nothing}
    T_in::Union{Float64,Nothing}
    fluid::Union{Fluid,Nothing}
    name::Union{String,Nothing}
    m_coolant::Union{Float64,Nothing}
    c_out::Union{Float64,Nothing}
end

function BreedingBlanket(;
    c_in=nothing, Q=nothing, TBR=nothing, T_out=nothing, T_in=nothing,
    fluid=nothing, name=nothing, m_coolant=nothing
)
    BreedingBlanket(c_in, Q, TBR, T_out, T_in, fluid, name, m_coolant, nothing)
end

"""Compute coolant mass flow rate from heat duty and temperature rise."""
function get_flowrate!(bb::BreedingBlanket)
    bb.m_coolant = bb.Q / ((bb.T_out - bb.T_in) * bb.fluid.cp)
end

"""
    get_cout!(bb; print_var=false)

Compute outlet tritium concentration from TBR, heat duty, and flow rate.
"""
function get_cout!(bb::BreedingBlanket; print_var::Bool=false)
    bb.m_coolant === nothing && get_flowrate!(bb)

    neutrons     = bb.Q / (REACTION_ENERGY * EV_TO_J)   # neutrons/s
    tritium_gen  = bb.TBR * neutrons / N_A               # mol/s

    if print_var
        println("neutrons/s: ", neutrons)
        println("T generated [mol/s]: ", tritium_gen)
    end

    vol_flow = bb.m_coolant / bb.fluid.rho  # volumetric flow [m³/s]
    if bb.fluid.MS
        bb.c_out = tritium_gen / 2 / vol_flow + bb.c_in   # T₂ in molten salt
    else
        bb.c_out = tritium_gen     / vol_flow + bb.c_in   # T  in liquid metal
    end
end

end # module BreedingBlanketModule
