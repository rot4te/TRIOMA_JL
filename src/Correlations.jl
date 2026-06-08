"""
    Correlations

Dimensionless-number and heat/mass-transfer correlations used across TRIOMA.
"""
module Correlations

export Nu_SiederTate, Nu_Gnielinsky, Nu_DittusBoelter, f_Pethukov, get_h_from_Nu,
       f_Haaland, Schmidt, Sherwood, Sherwood_HT_analogy, get_k_from_Sh, Re, Pr,
       Sherwood_bubbles, get_length_HX, get_deltaTML

"""Nusselt number — Sieder-Tate correlation."""
function Nu_SiederTate(Re, Pr, mu_w, mu_c)
    return 0.027 * Re^(4/5) * Pr^0.3 * (mu_c / mu_w)^0.14
end

"""Nusselt number — Gnielinski correlation."""
function Nu_Gnielinsky(Re, Pr, f)
    num = (f / 8) * (Re - 1000) * Pr
    den = 1 + 12.7 * (f / 8)^0.5 * (Pr^(2/3) - 1)
    return num / den
end

"""Nusselt number — Dittus-Boelter correlation."""
function Nu_DittusBoelter(Re, Pr)
    return 0.023 * Re^(4/5) * Pr^0.4
end

"""Friction factor — Pethukov correlation."""
function f_Pethukov(Re, Pr)
    return (0.79 * log(Re) - 1.64)^(-2)
end

"""Heat transfer coefficient from Nusselt number."""
function get_h_from_Nu(Nu, k, D)
    return Nu * k / D
end

"""Friction factor — Haaland correlation."""
function f_Haaland(Re, e_D)
    return (-1.8 * log10((e_D / 3.7)^1.11 + 6.9 / Re))^(-2)
end

"""Schmidt number."""
function Schmidt(D, mu, rho)
    return mu / (rho * D)
end

"""Sherwood number (Getthem paper correlation)."""
function Sherwood(Sc, Re)
    return 0.0096 * Re^0.913 * Sc^0.346
end

"""Sherwood number via heat-transfer analogy."""
function Sherwood_HT_analogy(Re, Sc)
    return 0.023 * Re^0.8 * Sc^0.4
end

"""Mass transfer coefficient from Sherwood number."""
function get_k_from_Sh(Sh, L, D)
    return Sh * D / L
end

"""Reynolds number."""
function Re(rho, u, L, mu)
    return rho * u * L / mu
end

"""Prandtl number."""
function Pr(c_p, mu, k)
    return mu * c_p / k
end

"""Sherwood number for bubbles (Humrickhouse MS report)."""
function Sherwood_bubbles(Sc, Re)
    return 0.089 * Re^0.69 * Sc^0.33
end

"""Length of a heat exchanger from LMTD, diameter, overall HTC, and heat duty."""
function get_length_HX(deltaTML, d_hyd, U, Q)
    return Q / (U * π * d_hyd * deltaTML)
end

"""Log-mean temperature difference for a counter-flow heat exchanger."""
function get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold)
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    return (delta_T1 - delta_T2) / log(delta_T1 / delta_T2)
end

end # module Correlations
