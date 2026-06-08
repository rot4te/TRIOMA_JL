"""
    ExtractorFunctions

Functions for calculating mass-transfer in packed-column Gas-Liquid Contactors (GLC).
Implements NTU integrals for both liquid-metal (Sievert) and molten-salt (Henry) systems.
"""
module ExtractorFunctions

using QuadGK
using Optim

export calculate_gas_velocity,
       extractor_lm, length_extractor_lm, NTU_lm, get_c_out_GLC_lm,
       extractor_ms, length_extractor_ms, NTU_ms, get_c_out_GLC_ms,
       pack_corr, corr_packed

const R_CONST = 8.314  # J / mol / K

"""Gas superficial velocity in the packed column [m/s]."""
function calculate_gas_velocity(G_gas, p_t, T, R)
    area  = π * R^2
    p_atm = 101325.0
    return G_gas / area / p_t * p_atm * T / 288.15
end

# ── Liquid-metal (Sievert) functions ──────────────────────────────────────────

"""
Number of transfer units (NTU) for liquid-metal systems.
Integrates the driving-force equation between c_out and c_in.
Returns the scalar integral value.
"""
function NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in; c_max=0.0)
    area  = π * R^2
    u_l   = G_l / area
    u_g   = calculate_gas_velocity(G_gas, p_t, T, R)
    R_g   = 2 * u_g / u_l
    c_in  = pl_in^0.5 * K_S
    c_out = pl_out^0.5 * K_S
    c_in_gas = pg_in / R_CONST / T

    toint(c) = 1 / (c - K_S * ((R_CONST * T / R_g) * (c - c_out + c_in_gas * R_g))^0.5)

    c_g_max   = pl_in / R_CONST / T
    c_out_max = max(
        c_in - R_g * (c_g_max - pg_in / R_CONST / T),
        c_max,
        pg_in^0.5 * K_S,
    )

    result, _ = quadgk(toint, c_out, c_in; atol=1e-10)
    if result < 0
        @warn "Negative NTU integral for LM: $(result)"
    end
    return result
end

"""Liquid load B_l [m/h] and k_la [1/s] from column height Z (liquid-metal)."""
function extractor_lm(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in)
    area  = π * R^2
    B_l   = G_l / area * 3600
    u_l   = G_l / area
    integ = NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in)
    kla_c = u_l / Z * integ
    return B_l, kla_c
end

"""Column height Z required for given kla (liquid-metal)."""
function length_extractor_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, kla; c_max=0.0)
    area  = π * R^2
    u_l   = G_l / area
    integ = NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in; c_max=c_max)
    return u_l / kla * integ
end

"""
    get_c_out_GLC_lm(Z, R, G_l, G_gas, pl_in, T, p_t, K_S, pg_in, kla)

Solve for outlet concentration (liquid-metal) given column height and kla.
Returns `(c_out, efficiency)`.
"""
function get_c_out_GLC_lm(Z, R, G_l, G_gas, pl_in, T, p_t, K_S, pg_in, kla)
    u_l   = G_l / (π * R^2)
    c_in  = pl_in^0.5 * K_S
    area  = π * R^2
    u_g   = calculate_gas_velocity(G_gas, p_t, T, R)

    c_out_max_reaction = c_in - kla * Z / u_l * (c_in - pg_in^0.5 * K_S)
    c_g_max = pl_in / R_CONST / T
    c_out_max = max(
        c_in - u_g / u_l * 2 * (c_g_max - pg_in / R_CONST / T),
        pg_in^0.5 * K_S,
        c_out_max_reaction,
    )

    function length_residual(c_out_v)
        c_out_s = c_out_v[1]
        pl_out_2 = c_out_s^2 / K_S^2
        z_guess = length_extractor_lm(R, G_l, G_gas, pl_in, pl_out_2, T, p_t, K_S, pg_in, kla; c_max=c_out_max)
        return abs(Z - z_guess)^2
    end

    x0   = [c_out_max + (c_in - c_out_max) / 2]
    lo   = [float(c_out_max)]
    hi   = [float(c_in)]
    res  = optimize(length_residual, lo, hi, x0, Fminbox(Powell()); options=Optim.Options(x_tol=1e-20, f_tol=1e-20, iterations=Int(1e8)))
    c_out = Optim.minimizer(res)[1]
    eff  = 1 - c_out / c_in

    L_cout = length_extractor_lm(R, G_l, G_gas, pl_in, c_out^2 / K_S^2, T, p_t, K_S, pg_in, kla)
    if abs(L_cout - Z) > 1e-3
        @warn "Guessed length ($(L_cout)) ≠ column height ($(Z)). Check your result."
    end
    if abs(c_out - c_out_max) < 1e-8
        println("The sweep gas saturated")
        if L_cout < Z
            println("Longer column would not increment the extraction efficiency")
        end
        eff = 1 - c_out / c_in
    end
    return c_out, eff
end

# ── Molten-salt (Henry) functions ─────────────────────────────────────────────

"""
Number of transfer units (NTU) for molten-salt systems.
"""
function NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in; c_max=0.0)
    area     = π * R^2
    u_l      = G_l / area
    c_in     = pl_in  * K_H
    c_out    = pl_out * K_H
    u_g      = calculate_gas_velocity(G_gas, p_t, T, R)
    c_in_gas = pg_in / R_CONST / T
    R_g      = u_g / u_l

    toint(c) = 1 / (c - K_H * (u_l / u_g * R_CONST * T) * (c - c_out + c_in_gas * u_g / u_l))

    # Python source uses fixed_quad for the final result; use quadgk here.
    result, _ = quadgk(toint, c_out, c_in; atol=1e-10)
    return result
end

"""Liquid load B_l [m/h] and k_la [1/s] from column height Z (molten-salt)."""
function extractor_ms(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in)
    area  = π * R^2
    B_l   = G_l / area * 3600
    u_l   = G_l / area
    integ = NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in)
    kla_c = u_l / Z * integ
    return B_l, kla_c
end

"""Column height Z required for given kla (molten-salt)."""
function length_extractor_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, kla; c_max=0.0)
    area  = π * R^2
    u_l   = G_l / area
    integ = NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in; c_max=c_max)
    return u_l / kla * integ
end

"""
    get_c_out_GLC_ms(Z, R, G_l, G_gas, pl_in, T, p_t, K_H, pg_in, kla)

Solve for outlet concentration (molten-salt) given column height and kla.
Returns `(c_out, efficiency)`.
"""
function get_c_out_GLC_ms(Z, R, G_l, G_gas, pl_in, T, p_t, K_H, pg_in, kla)
    u_l  = G_l / (π * R^2)
    c_in = pl_in * K_H
    u_g  = calculate_gas_velocity(G_gas, p_t, T, R)

    c_out_max_reaction = c_in - kla * Z / u_l * (c_in - pg_in * K_H)
    c_g_max  = pl_in / R_CONST / T
    c_out_max = max(
        c_in - u_g / u_l * (c_g_max - pg_in / R_CONST / T),
        pg_in * K_H,
        c_out_max_reaction,
    )

    function length_residual(c_out_v)
        c_out_s  = c_out_v[1]
        pl_out_2 = c_out_s / K_H
        z_guess  = length_extractor_ms(R, G_l, G_gas, pl_in, pl_out_2, T, p_t, K_H, pg_in, kla; c_max=c_out_max)
        return abs(Z - z_guess)^2
    end

    x0  = [c_out_max + (c_in - c_out_max) / 2]
    lo  = [float(c_out_max)]
    hi  = [float(c_in)]
    res = optimize(length_residual, lo, hi, x0, Fminbox(Powell()); options=Optim.Options(x_tol=1e-20, f_tol=1e-20, iterations=Int(1e8)))
    c_out = Optim.minimizer(res)[1]
    eff   = 1 - c_out / c_in

    L_cout = length_extractor_ms(R, G_l, G_gas, pl_in, c_out / K_H, T, p_t, K_H, pg_in, kla)
    if abs(L_cout - Z) > 1e-3
        @warn "Guessed length ($(L_cout)) ≠ column height ($(Z)). Check your result."
    end
    if abs(c_out - c_out_max) < 1e-8
        println("The sweep gas saturated")
        if L_cout < Z
            println("Longer column would not increment the extraction efficiency")
        end
        eff = 1 - c_out / c_in
    end
    return c_out, eff
end

# ── Packed-column correlations ────────────────────────────────────────────────

"""
Liquid-film mass transfer coefficient for a packed column (Onda correlation).
"""
function pack_corr(a, d, D, eta, v)
    k_l = 0.0051 * (v / eta / a)^(2/3) * (D / eta)^0.5 * (a * d)^0.4 * (1 / eta / 9.81)^(-1/3)
    return k_l
end

"""
Sherwood-based k_l for packed column with Raschig rings.

Warning: verification of this correlation is noted as incomplete in the
original Python source.
"""
function corr_packed(Re, Sc, d, rho_L, mu_L, L, D)
    using ..Correlations: get_k_from_Sh
    beta = 0.32  # Raschig rings (0.25 for some references)
    g    = 9.81
    Sh   = beta * Re^0.59 * Sc^0.5 * (d^3 * g * rho_L^2 / mu_L^2)^0.17
    return get_k_from_Sh(Sh, L, D)
end

end # module ExtractorFunctions
