"""
    Correlations

Dimensionless-number and heat/mass-transfer correlations used across TRIOMA.
"""
module Correlations

export Nu_SiederTate, Nu_Gnielinsky, Nu_DittusBoelter, f_Pethukov, get_h_from_Nu,
       f_Haaland, Schmidt, Sherwood, Sherwood_HT_analogy, get_k_from_Sh, Re, Pr,
       Sherwood_bubbles, get_length_HX, get_deltaTML, hydraulic_diameter,
       Nu_Iev_tube, Nu_Asma_tube, Nu_Si_tube,
       Nu_Yang_lam, Nu_Yang1_tube, Nu_Yang2_tube,
       fD_Iev_tube, fD_Asma_tube, fD_Si_tube, fD_Gao_tube, fD_Yang_tube

"""
    hydraulic_diameter(area, wetted_perimeter)

Hydraulic diameter of a flow cross-section, ``D_h = 4 A / P``, where `area` is
the wetted (flow) cross-sectional area and `wetted_perimeter` is the wetted
perimeter. This is the general definition valid for any cross-section; for a
full circular pipe it reduces to the inner diameter.
"""
hydraulic_diameter(area::Float64, wetted_perimeter::Float64)::Float64 =
    4 * area / wetted_perimeter

###########################################################################
# Twisted-tube (TET) correlations — tube-side (internal flow)
# Ported from HeatExchange/src/TTHX-correlations.jl.
# D_max = major inner diameter (2a), D_min = minor inner diameter (2b),
# D_h   = hydraulic diameter (4A/P),  s = twist pitch [m] (axial length
# per full 360° rotation).  ϕ (viscosity wall-correction) is omitted (= 1).
###########################################################################

"""
    Nu_Iev_tube(Re, D_max, s)

Tube-side Nusselt number for a twisted elliptical tube — Ievlev correlation.
No Prandtl dependence (correlation developed for a specific fluid regime).

    Nu = 0.019 Re^0.8 (1 + 0.547 (s/D_max)^−0.83)

Valid: `6000 ≤ Re ≤ 100 000`, `s/D_max ≥ 6.2`.
"""
Nu_Iev_tube(Re::Float64, D_max::Float64, s::Float64)::Float64 =
    0.019 * Re^0.8 * (1 + 0.547 * (s / D_max)^(-0.83))

"""
    Nu_Asma_tube(Re, Pr, D_max, s)

Tube-side Nusselt number for a twisted elliptical tube — Asmantas Modified
correlation (Hughes 2017, eq. A.8a), which replaces the original gas-phase
`(Tw/Tf)^n` viscosity correction with `(µw/µ)^−0.14`. Here ϕ = 1 (no
wall-correction; use when bulk ≈ wall viscosity or as a first estimate).

    Nu = 0.021 Re^0.8 Pr^0.4 (1 + 2.1 (s/D_max)^−0.91)

Valid: `7000 ≤ Re ≤ 200 000`, `6.2 ≤ s/D_max ≤ 12.2`.
"""
Nu_Asma_tube(Re::Float64, Pr::Float64, D_max::Float64, s::Float64)::Float64 =
    0.021 * Re^0.8 * Pr^0.4 * (1 + 2.1 * (s / D_max)^(-0.91))

"""
    Nu_Si_tube(Re, Pr, D_h, D_max, s)

Tube-side Nusselt number for a twisted elliptical tube — Simonin correlation.

    Nu = 0.396 Re^0.544 (s/D_h)^0.161 (s/D_max)^−0.519 Pr^0.33

Valid: `1000 ≤ Re ≤ 17 000`, `6.86 ≤ s/D_max ≤ 11.9`, `0.144 ≤ s ≤ 0.250 m`.
"""
Nu_Si_tube(Re::Float64, Pr::Float64, D_h::Float64, D_max::Float64, s::Float64)::Float64 =
    0.396 * Re^0.544 * (s / D_h)^0.161 * (s / D_max)^(-0.519) * Pr^0.33

"""
    Nu_Yang1_tube(Re, Pr, D_h, D_max, D_min, s)

Tube-side Nusselt number for a twisted elliptical tube — Yang correlation 1.

    Nu = 0.034 Re^0.784 Pr^0.333 (D_min/D_max)^−0.590 (s/D_h)^−0.165

Valid: `5000 ≤ Re ≤ 20 000`, `0.2 ≤ s ≤ 0.4 m`,
`0.0198 ≤ D_max ≤ 0.0218 m`, `0.0058 ≤ D_min ≤ 0.0094 m`.
"""
Nu_Yang1_tube(Re::Float64, Pr::Float64, D_h::Float64,
              D_max::Float64, D_min::Float64, s::Float64)::Float64 =
    0.034 * Re^0.784 * Pr^0.333 * (D_min / D_max)^(-0.590) * (s / D_h)^(-0.165)

"""
    Nu_Yang2_tube(Re, Pr, D_h, D_max, D_min, s)

Tube-side Nusselt number for a twisted elliptical tube — Yang correlation 2.

    Nu = 1.50618 Re^0.51825 Pr^−1.2446 (D_max/D_min)^1.12252 (s/D_h)^−0.32367

Note: the negative Pr exponent is as published; Pr appears in the denominator.
Valid: `7900 ≤ Re ≤ 26 500`, `0.160 ≤ s ≤ 0.250 m`, `0.024 ≤ D_max ≤ 0.026 m`.
"""
Nu_Yang2_tube(Re::Float64, Pr::Float64, D_h::Float64,
              D_max::Float64, D_min::Float64, s::Float64)::Float64 =
    1.50618 * Re^0.51825 * Pr^(-1.2446) * (D_max / D_min)^1.12252 * (s / D_h)^(-0.32367)

"""
    Nu_Yang_lam(Re, Pr, D_max, D_min, s)

Tube-side Nusselt number for a twisted elliptical tube in laminar flow — Yang
Laminar correlation (Hughes 2017, eq. A.10a). Hughes' Table 4.3 found this the
best-performing tube-side correlation across forced and natural circulation
(MARE 7.7 %).

    Nu = 3.66 + 0.512 Re^0.477 Pr^0.975 (1 − D_min/D_max)^1.532 (s/D_min)^−0.609

Note: uses `D_min` (minor diameter = 2b) in the twist-pitch ratio, not `D_h`.
Valid: `100 ≤ Re ≤ 500`, `10 ≤ s/D_max ≤ 17`.
"""
Nu_Yang_lam(Re::Float64, Pr::Float64, D_max::Float64, D_min::Float64, s::Float64)::Float64 =
    3.66 + 0.512 * Re^0.477 * Pr^0.975 * (1 - D_min / D_max)^1.532 * (s / D_min)^(-0.609)

###########################################################################
# Twisted-tube friction factors — tube-side
###########################################################################

"""
    fD_Iev_tube(Re, D_max, s)

Tube-side Darcy friction factor for a twisted elliptical tube — Ievlev correlation.

    f = 0.316 (1 + 3.27 (s/D_max)^−0.87) Re^−0.25

Valid: `6000 ≤ Re ≤ 100 000`, `6.2 ≤ s/D_max ≤ 16.7`.
"""
fD_Iev_tube(Re::Float64, D_max::Float64, s::Float64)::Float64 =
    0.316 * (1 + 3.27 * (s / D_max)^(-0.87)) * Re^(-0.25)

"""
    fD_Asma_tube(Re, D_max, s)

Tube-side Darcy friction factor for a twisted elliptical tube — Asmaa correlation.

    f = 0.82 (s/D_max)^−0.63 Re^−0.18

Valid: `7000 ≤ Re ≤ 200 000`, `6.2 ≤ s/D_max ≤ 12.2`.
"""
fD_Asma_tube(Re::Float64, D_max::Float64, s::Float64)::Float64 =
    0.82 * (s / D_max)^(-0.63) * Re^(-0.18)

"""
    fD_Si_tube(Re, D_max, s)

Tube-side Darcy friction factor for a twisted elliptical tube — Simonin correlation
(Hughes 2017, eq. A.16).

    f = 10^(a1 + a2 log₁₀(Re) + a3 log₁₀(Re)²)

where a1–a3 are polynomials in s/D_max (eqs. A.16b–d). Hughes' Appendix C Python
code uses `log10` throughout; the `log` in eq. A.16 means base-10.

Valid: `1000 ≤ Re ≤ 17 000`.
"""
function fD_Si_tube(Re::Float64, D_max::Float64, s::Float64)::Float64
    r  = s / D_max
    a1 = -19.70 + 4.90*r - 0.22*r^2
    a2 =  10.52 - 2.66*r + 0.12*r^2
    a3 =  -1.47 + 0.36*r - 0.016*r^2
    return 10^(a1 + a2*log10(Re) + a3*log10(Re)^2)
end

"""
    fD_Gao_tube(Re, D_h, D_max, D_min, s)

Tube-side Darcy friction factor for a twisted elliptical tube — Gao correlation.

    f = 4.572 Re^−0.521 (D_min/D_max)^−0.334 (s/D_h)^−0.082

Valid: `5000 ≤ Re ≤ 20 000`.
"""
fD_Gao_tube(Re::Float64, D_h::Float64, D_max::Float64, D_min::Float64, s::Float64)::Float64 =
    4.572 * Re^(-0.521) * (D_min / D_max)^(-0.334) * (s / D_h)^(-0.082)

"""
    fD_Yang_tube(Re, Pr, D_h, D_max, D_min, s)

Tube-side Darcy friction factor for a twisted elliptical tube — Yang correlation
(Hughes 2017, eq. A.17a). Note: the HeatExchange source omitted `s` from the
parameter list; corrected here per the dissertation.

    f = 0.71497 Re^0.07777 Pr^−1.03974 (D_max/D_min)^−0.076212 (s/D_h)^−0.33393

Valid: `7900 ≤ Re ≤ 26 500`, `0.160 ≤ s ≤ 0.250 m`, `0.024 ≤ D_max ≤ 0.026 m`.
"""
fD_Yang_tube(Re::Float64, Pr::Float64, D_h::Float64,
             D_max::Float64, D_min::Float64, s::Float64)::Float64 =
    0.71497 * Re^0.07777 * Pr^(-1.03974) * (D_max / D_min)^(-0.076212) * (s / D_h)^(-0.33393)

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
