"""
BreedingBlanket.jl

Tritium-producing breeding blanket component.
Computes the tritium outlet concentration based on fusion power, TBR,
and coolant flow rate.
"""

module BreedingBlanketModule

using ..TriomaTypes
using ..PipeSubclasses: Fluid, FluidMaterial

export BreedingBlanket, get_flowrate!, get_cout!, connect_to_component!

# Physical constants (matching scipy.constants values used in Python source)
using AtomicAndPhysicalConstants: J_PER_EV, AVOGADRO

const eV_TO_J = J_PER_EV   # J per eV
const N_A     = AVOGADRO    # mol⁻¹
const REACTION_ENERGY_EV = 17.6e6  # eV  (DT fusion: 17.6 MeV)

# ---------------------------------------------------------------------------
# BreedingBlanket struct
# ---------------------------------------------------------------------------

"""
Breeding blanket component.

Produces tritium from fusion neutrons; computes the tritium outlet
concentration `c_out` [mol/m³] given fusion power `Q`, tritium breeding
ratio `TBR`, and coolant properties.

Fields:
- `c_in`      : inlet tritium concentration [mol/m³]
- `c_out`     : outlet tritium concentration [mol/m³] (set by `get_cout!`)
- `Q`         : fusion thermal power [W]
- `TBR`       : tritium breeding ratio (dimensionless)
- `T_out`     : coolant outlet temperature [K]
- `T_in`      : coolant inlet temperature [K]
- `fluid`     : `Fluid` object representing the coolant
- `name`      : optional label string
- `m_coolant` : coolant mass flow rate [kg/s] (computed by `get_flowrate!`)
"""
mutable struct BreedingBlanket
    c_in     ::Union{Float64, Nothing}
    c_out    ::Union{Float64, Nothing}
    Q        ::Union{Float64, Nothing}
    TBR      ::Union{Float64, Nothing}
    T_out    ::Union{Float64, Nothing}
    T_in     ::Union{Float64, Nothing}
    fluid    ::Union{Fluid, FluidMaterial, Nothing}
    name     ::Union{String, Nothing}
    m_coolant::Union{Float64, Nothing}
end

BreedingBlanket(;
    c_in=nothing, c_out=nothing, Q=nothing, TBR=nothing,
    T_out=nothing, T_in=nothing, fluid=nothing,
    name=nothing, m_coolant=nothing) =
    BreedingBlanket(c_in, c_out, Q, TBR, T_out, T_in, fluid, name, m_coolant)

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

"""
    get_flowrate!(bb::BreedingBlanket)

Compute the coolant mass flow rate [kg/s] from `Q`, `T_out - T_in`, and
`fluid.cp`. Stores the result in `bb.m_coolant`.
"""
function get_flowrate!(bb::BreedingBlanket)
    bb.m_coolant = bb.Q / ((bb.T_out - bb.T_in) * bb.fluid.cp)
end

"""
    get_cout!(bb::BreedingBlanket; print_var=false)

Compute the tritium outlet concentration [mol/m³] and store it in `bb.c_out`.

Physics:
- Number of DT fusion reactions per second: `neutrons = Q / (17.6 MeV × eV_to_J)`
- Tritium produced per second: `tritium_gen = TBR × neutrons / Nₐ`  [mol/s]
- For molten-salt (MS=true):  tritium exists as T₂, so divide by 2.
  `c_out = tritium_gen / 2 / (m_coolant / ρ) + c_in`
- For liquid metal (MS=false): tritium is atomic T.
  `c_out = tritium_gen     / (m_coolant / ρ) + c_in`
"""
function get_cout!(bb::BreedingBlanket; print_var::Bool=false)
    bb.m_coolant === nothing && get_flowrate!(bb)

    neutrons    = bb.Q / (REACTION_ENERGY_EV * eV_TO_J)
    tritium_gen = bb.TBR * neutrons / N_A   # mol/s

    if print_var
        println("neutrons    = ", neutrons)
        println("tritium_gen = ", tritium_gen, " mol/s")
    end

    vol_flow = bb.m_coolant / bb.fluid.rho   # m³/s

    if bb.fluid.MS
        bb.c_out = tritium_gen / 2 / vol_flow + bb.c_in
    else
        bb.c_out = tritium_gen     / vol_flow + bb.c_in
    end
end

"""
    connect_to_component!(bb::BreedingBlanket, next)

Set the inlet concentration of `next` to the outlet concentration of `bb`.
`next` can be any TRIOMA component that has a `c_in` field.
"""
function connect_to_component!(bb::BreedingBlanket, next)
    next === nothing && throw(ArgumentError("next component cannot be nothing"))
    update_attribute!(next, "c_in", bb.c_out)
end

end # module
