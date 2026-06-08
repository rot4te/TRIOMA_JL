"""
    MoltenSalts

Dimensionless parameters and regime identification for molten-salt breeders.
"""
module MoltenSalts

export W_ms, H_ms, get_regime_ms

"""
W parameter (surface/diffusion ratio) for molten salts.

W << 1 → surface limited; W >> 1 → diffusion limited.
"""
function W_ms(k_d, D, thick, K_S, c0, k_H)
    return 2 * k_d * thick * (c0 / k_H)^0.5 / (D * K_S)
end

"""
H dimensionless number for molten salts.

H >> 1 and H/W >> 1 → mass-transport limited.
"""
function H_ms(k_t, k_H, k_d)
    return k_d / (k_t * k_H)
end

"""
    get_regime_ms(; k_d, D, thick, K_S, c0, k_t, k_H, print_var=false)

Return a string describing the dominant transport regime for a molten-salt system.
"""
function get_regime_ms(; k_d, D, thick, K_S, c0, k_t, k_H, print_var::Bool=false)
    W = W_ms(k_d, D, thick, K_S, c0, k_H)
    H = H_ms(k_t, k_H, k_d)
    result = if H > 10 && H / W > 10
        "Mass transport limited"
    elseif H < 0.1 && W < 0.1
        "Surface limited"
    elseif W > 10 && H / W < 0.1
        "Diffusion Limited"
    else
        "Mixed regime"
    end
    if print_var
        println("H is equal to ", H)
        println("W is equal to ", W)
        println("H/W is equal to ", H / W)
        println(result)
    end
    return result
end

end # module MoltenSalts


"""
    LiquidMetals

Dimensionless parameters and regime identification for liquid-metal breeders.
"""
module LiquidMetals

export W_lm, partition_param_lm, get_regime_lm

"""
W parameter (surface/diffusion ratio) for liquid metals.

W << 1 → surface limited.
"""
function W_lm(k_r, D, thick, K_S, c0, K_S_L)
    return k_r / D * K_S * thick * c0 / K_S_L
end

"""
Partition parameter (diffusion vs mass-transfer) for liquid metals.
"""
function partition_param_lm(D, k_t, K_S_S, K_S_L, t)
    return D / k_t * K_S_S / (K_S_L * t)
end

"""
    get_regime_lm(; D, k_t, K_S_S, K_S_L, k_r, thick, c0, print_var=false)

Return a string describing the dominant transport regime for a liquid-metal system.
"""
function get_regime_lm(; D, k_t, K_S_S, K_S_L, k_r, thick, c0, print_var::Bool=false)
    W   = W_lm(k_r, D, thick, K_S_S, c0, K_S_L)
    pp  = partition_param_lm(D, k_t, K_S_S, K_S_L, thick)
    result = if pp > 10 && W > 10
        "Mass transport limited"
    elseif pp < 0.1 && W < 0.1
        "Surface limited"
    elseif W > 10 && pp < 0.1
        "Diffusion Limited"
    elseif pp > 10 && W < 0.1
        "Transport and surface limited regime"
    else
        "Mixed regime"
    end
    if print_var
        println("H is equal to ", pp * W)
        println("W is equal to ", W)
        println("H/W is equal to ", pp)
        println(result)
    end
    return result
end

end # module LiquidMetals
