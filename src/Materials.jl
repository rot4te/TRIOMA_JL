"""
    Materials

Pre-defined material property functions for common tritium-breeding fluids and
structural metals. Each function returns a `FluidMaterial` or `SolidMaterial`
evaluated at temperature `T` [K].

Note: several LiPb and Sodium properties are currently placeholder values
(marked with # TODO in the original Python source).
"""
module Materials

using ..PipeSubclasses: FluidMaterial, SolidMaterial
using AtomicAndPhysicalConstants: BOLTZMANN_k, J_PER_EV, AVOGADRO

export Flibe, Sodium, LiPb, Steel

const R_CONST = BOLTZMANN_k * J_PER_EV * AVOGADRO   # J / mol / K

# ── Flibe ─────────────────────────────────────────────────────────────────────

function Flibe(T::Real)::FluidMaterial
    density(T)    = 2413 - 0.488 * T
    viscosity(T)  = 1.16e-4 * exp(3755 / T)
    H_diff(T)     = 9.3e-7 * exp(-42e3 / (R_CONST * T))
    k_thermal()   = 1.1     # W/m/K
    cp_val()      = 2386.0  # J/kg/K
    k_H(T)        = 4.54e-4 # Henry's constant mol/m³/Pa (constant fit)

    return FluidMaterial(
        T         = Float64(T),
        D         = H_diff(T),
        Solubility = k_H(T),
        MS        = true,
        mu        = viscosity(T),
        rho       = density(T),
        k         = k_thermal(),
        cp        = cp_val(),
    )
end

# ── Sodium ────────────────────────────────────────────────────────────────────

function Sodium(T::Real)::FluidMaterial
    density(T)   = 219 + 275.32 * (1 - T / 2504.7) + 511.58 * (1 - T / 2503.7)^0.5
    viscosity(T) = exp(-6.4406 - 0.3958 * log(T) + 556.835 / T)
    H_diff(T)    = 2e-5 * exp(-49053 / (R_CONST * T))
    k_S(T)       = 10^(0.86 - 122 / T)  # Sievert constant mol/m³/Pa

    # Note: k and cp are placeholder (TODO in original source)
    return FluidMaterial(
        T         = Float64(T),
        D         = H_diff(T),
        Solubility = k_S(T),
        MS        = false,
        mu        = viscosity(T),
        rho       = density(T),
        k         = 1.0,   # TODO
        cp        = 1.0,   # TODO
    )
end

# ── LiPb ─────────────────────────────────────────────────────────────────────

function LiPb(T::Real)::FluidMaterial
    density(T)   = 9659.8  # TODO
    viscosity(T) = 1.0     # TODO
    H_diff(T)    = 1.0     # TODO
    MM_LiPb      = 180.0   # TODO
    k_S(T)       = 4.7e-7 * exp(-9e3 / (R_CONST * T) * density(T) / MM_LiPb)  # TODO

    return FluidMaterial(
        T         = Float64(T),
        D         = H_diff(T),
        Solubility = k_S(T),
        MS        = false,
        mu        = viscosity(T),
        rho       = density(T),
        k         = 1.0,   # TODO
        cp        = 1.0,   # TODO
    )
end

# ── Steel ─────────────────────────────────────────────────────────────────────

function Steel(T::Real)::SolidMaterial
    H_diff(T) = 5.81e-7 * exp(-66.3e3 / (R_CONST * T))
    K_S_val   = 1.0  # placeholder

    return SolidMaterial(
        T   = Float64(T),
        D   = H_diff(T),
        K_S = K_S_val,
    )
end

end # module Materials
