"""
    Correlations

Dimensionless-number and heat/mass-transfer correlations used across TRIOMA.
"""
module Correlations

export Nu_SiederTate, Nu_Gnielinsky, Nu_DittusBoelter, f_Pethukov, get_h_from_Nu,
       f_Haaland, Schmidt, Sherwood, Sherwood_HT_analogy, get_k_from_Sh, Re, Pr,
       Sherwood_bubbles, get_length_HX, get_deltaTML

"""Nusselt number — Sieder-Tate correlation."""
function Nu_SiederTate(Re::Float64, Pr::Float64, mu_w::Float64, mu_c::Float64)::Float64
    return 0.027 * Re^(4/5) * Pr^0.3 * (mu_c / mu_w)^0.14
end

"""Nusselt number — Gnielinski correlation."""
function Nu_Gnielinsky(Re::Float64, Pr::Float64, f::Float64)::Float64
    num = (f / 8) * (Re - 1000) * Pr
    den = 1 + 12.7 * (f / 8)^0.5 * (Pr^(2/3) - 1)
    return num / den
end

"""Nusselt number — Dittus-Boelter correlation."""
Nu_DittusBoelter(Re::Float64, Pr::Float64)::Float64 = 0.023 * Re^(4/5) * Pr^0.4

"""Friction factor — Pethukov correlation."""
f_Pethukov(Re::Float64, ::Float64)::Float64 = (0.79 * log(Re) - 1.64)^(-2)

"""Heat transfer coefficient from Nusselt number."""
get_h_from_Nu(Nu::Float64, k::Float64, D::Float64)::Float64 = Nu * k / D

"""Friction factor — Haaland correlation."""
function f_Haaland(Re::Float64, e_D::Float64)::Float64
    return (-1.8 * log10((e_D / 3.7)^1.11 + 6.9 / Re))^(-2)
end

"""Schmidt number."""
Schmidt(D::Float64, mu::Float64, rho::Float64)::Float64 = mu / (rho * D)

"""Sherwood number (Getthem paper correlation)."""
Sherwood(Sc::Float64, Re::Float64)::Float64 = 0.0096 * Re^0.913 * Sc^0.346

"""Sherwood number via heat-transfer analogy."""
Sherwood_HT_analogy(Re::Float64, Sc::Float64)::Float64 = 0.023 * Re^0.8 * Sc^0.4

"""Mass transfer coefficient from Sherwood number."""
get_k_from_Sh(Sh::Float64, L::Float64, D::Float64)::Float64 = Sh * D / L

"""Reynolds number."""
Re(rho::Float64, u::Float64, L::Float64, mu::Float64)::Float64 = rho * u * L / mu

"""Prandtl number."""
Pr(c_p::Float64, mu::Float64, k::Float64)::Float64 = mu * c_p / k

"""Sherwood number for bubbles (Humrickhouse MS report)."""
Sherwood_bubbles(Sc::Float64, Re::Float64)::Float64 = 0.089 * Re^0.69 * Sc^0.33

"""Length of a heat exchanger from LMTD, diameter, overall HTC, and heat duty."""
function get_length_HX(deltaTML::Float64, d_hyd::Float64, U::Float64, Q::Float64)::Float64
    return Q / (U * π * d_hyd * deltaTML)
end

"""Log-mean temperature difference for a counter-flow heat exchanger."""
function get_deltaTML(T_in_hot::Float64, T_out_hot::Float64,
                      T_in_cold::Float64, T_out_cold::Float64)::Float64
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    return (delta_T1 - delta_T2) / log(delta_T1 / delta_T2)
end

end # module Correlations
