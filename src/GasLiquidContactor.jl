"""
    GasLiquidContactor

Types and methods for Gas-Liquid Contactor (packed-column) components.
"""
module GasLiquidContactor

using ..TriomaTypes: update_attribute!
using ..PipeSubclasses: Fluid, Membrane
using ..ExtractorFunctions

export GLC_Gas, GLC,
       get_c_out!, get_kla_from_cout!, get_z_from_eff

# ── GLC_Gas ───────────────────────────────────────────────────────────────────

"""
    GLC_Gas

Sweep gas stream for a Gas-Liquid Contactor.

Fields
- `G_gas`  : volumetric flow rate at normal conditions [m³/s]
- `pg_in`  : inlet tritium partial pressure [Pa] (default 0)
- `pg_out` : outlet tritium partial pressure [Pa] (default 0)
- `p_tot`  : total column pressure [Pa] (default 100 000)
"""
mutable struct GLC_Gas
    G_gas::Union{Float64,Nothing}
    pg_in::Float64
    pg_out::Float64
    p_tot::Float64
end
GLC_Gas(; G_gas=nothing, pg_in=0.0, pg_out=0.0, p_tot=100_000.0) =
    GLC_Gas(G_gas, pg_in, pg_out, p_tot)

# ── GLC ───────────────────────────────────────────────────────────────────────

"""
    GLC

Gas-Liquid Contactor (packed column) component.

Fields
- `H`       : column height [m]
- `R`       : column radius [m]
- `L`       : characteristic length for fluid flow [m]
- `c_in`    : inlet tritium concentration [mol/m³]
- `c_out`   : outlet tritium concentration [mol/m³]
- `eff`     : extraction efficiency
- `fluid`   : liquid-phase `Fluid` object
- `GLC_gas` : gas-phase `GLC_Gas` object
- `T`       : operating temperature [K]
- `G_L`     : liquid volumetric flow rate [m³/s]
- `kla`     : overall mass-transfer coefficient × packing area [1/s]
"""
mutable struct GLC
    H::Union{Float64,Nothing}
    R::Union{Float64,Nothing}
    L::Union{Float64,Nothing}
    c_in::Union{Float64,Nothing}
    c_out::Union{Float64,Nothing}
    eff::Union{Float64,Nothing}
    fluid::Union{Fluid,Nothing}
    GLC_gas::Union{GLC_Gas,Nothing}
    T::Union{Float64,Nothing}
    G_L::Union{Float64,Nothing}
    kla::Union{Float64,Nothing}
end

function GLC(;
    H=nothing, R=nothing, L=nothing, c_in=nothing, c_out=nothing, eff=nothing,
    fluid=nothing, GLC_gas=nothing, T=nothing, G_L=nothing, kla=nothing
)
    GLC(H, R, L, c_in, c_out, eff, fluid, GLC_gas, T, G_L, kla)
end

# ── Methods ───────────────────────────────────────────────────────────────────

"""
    get_c_out!(glc::GLC)

Compute `glc.c_out` and `glc.eff` from current state.
Dispatches on `glc.fluid.MS` (molten salt vs liquid metal).
"""
function get_c_out!(glc::GLC)
    if !glc.fluid.MS
        c_out, eff = ExtractorFunctions.get_c_out_GLC_lm(
            glc.H, glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in^2 / glc.fluid.Solubility^2,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility,
            glc.GLC_gas.pg_in, glc.kla
        )
    else
        c_out, eff = ExtractorFunctions.get_c_out_GLC_ms(
            glc.H, glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in / glc.fluid.Solubility,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility,
            glc.GLC_gas.pg_in, glc.kla
        )
    end
    glc.c_out = c_out
    glc.eff   = eff
end

"""
    get_kla_from_cout!(glc::GLC) -> (Bl, kla)

Compute `kla` from a known `c_out`. Sets `glc.kla` and `glc.Bl`.
"""
function get_kla_from_cout!(glc::GLC)
    if !glc.fluid.MS
        Bl, kla = ExtractorFunctions.extractor_lm(
            glc.H, glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in^2  / glc.fluid.Solubility^2,
            glc.c_out^2 / glc.fluid.Solubility^2,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility, glc.GLC_gas.pg_in
        )
    else
        Bl, kla = ExtractorFunctions.extractor_ms(
            glc.H, glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in  / glc.fluid.Solubility,
            glc.c_out / glc.fluid.Solubility,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility, glc.GLC_gas.pg_in
        )
    end
    glc.kla = kla
    return Bl, kla
end

"""
    get_z_from_eff(glc::GLC) -> Z

Return the column height needed to achieve the current `glc.eff`.
"""
function get_z_from_eff(glc::GLC)
    if !glc.fluid.MS
        return ExtractorFunctions.length_extractor_lm(
            glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in^2  / glc.fluid.Solubility^2,
            glc.c_out^2 / glc.fluid.Solubility^2,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility, glc.GLC_gas.pg_in, glc.kla
        )
    else
        return ExtractorFunctions.length_extractor_ms(
            glc.R, glc.G_L, glc.GLC_gas.G_gas,
            glc.c_in  / glc.fluid.Solubility,
            glc.c_out / glc.fluid.Solubility,
            glc.T, glc.GLC_gas.p_tot, glc.fluid.Solubility, glc.GLC_gas.pg_in, glc.kla
        )
    end
end

end # module GasLiquidContactor
