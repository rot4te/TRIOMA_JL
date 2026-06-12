"""
    PAVPipe

Permeation Against Vacuum (PAV) pipe component: the core building block of
TRIOMA's outer fuel cycle model.

Exports the `Component` mutable struct and all functions that operate on it,
covering:
- Hydraulics: friction factor, pressure drop, pumping power
- Transport-regime identification and dimensionless parameters
- Analytical and numerical extraction efficiency
- Permeation flux (all regimes, both MS and LM chemistries)
- Tritium inventory in fluid and membrane
- Heat-exchanger overall coefficient
"""
module PAVPipe

using ..TriomaCore
import ..TriomaCore: update_attribute!
using ..PipeSubclasses
using ..Correlations
using ..FusionCoolant
using LambertW
using QuadGK
using Optim

export Component,
       get_regime, get_pipe_flowrate, get_total_flowrate,
       define_component_volumes!, get_adimensionals!,
       use_analytical_efficiency!, analytical_efficiency!,
       get_efficiency!, get_flux!, outlet_c_comp!,
       get_global_HX_coeff!, friction_factor, get_pressure_drop!,
       get_pumping_power!, estimate_cost!, set_hydraulic_diameter!,
       get_solid_inventory!, get_fluid_inventory!, get_inventory!,
       analytical_solid_inventory!, analytical_fluid_inventory!,
       T_leak, update_T_prop!

# ---------------------------------------------------------------------------
# Component struct
# ---------------------------------------------------------------------------

"""
    Component

A single pipe-type component in the TRIOMA outer-fuel-cycle model. Covers
both permeation-against-vacuum (PAV) extractors and heat-exchanger pipes.

All physical inputs are set by the user (or copied from sub-objects); all
derived quantities are written back by the calculation functions.

## Required inputs
- `geometry`  : `Geometry` with inner diameter `D`, wall thickness `dw`, length `L`, and `n_pipes`
- `fluid`     : `Fluid` with tritium diffusivity `D`, solubility, velocity `U0`, and thermo-physical properties
- `membrane`  : `Membrane` with diffusivity `D`, Sievert constant `K_S`, surface rate `k_d`, wall thickness `dw`
- `c_in`      : tritium inlet concentration [mol/m³]

## Computed outputs (set by calculation methods)
- `c_out`     : outlet concentration [mol/m³]
- `eff`       : extraction efficiency (numerical)
- `eff_an`    : extraction efficiency (analytical)
- `J_perm`    : local permeation flux [mol/m²/s]  (negative = extraction from fluid)
- `H`, `W`    : dimensionless transport parameters
- `tau`       : dimensionless residence time
- `alpha`/`xi`: MS permeability parameters
- `zeta`      : LM permeability parameter
- `inv`       : total tritium inventory [mol]
- `delta_p`   : pressure drop [Pa]
- `U`         : overall heat-transfer coefficient [W/m²/K]
- `pumping_power` : pumping power [W]

## Other fields
- `name`      : optional label
- `loss`      : if `true`, this component counts as a loss pipe in circuit analysis
- `p_out`     : downstream tritium partial pressure [Pa] (default 1e-15, near vacuum)
"""
mutable struct Component
    geometry      ::Union{Geometry, Nothing}
    c_in          ::Union{Float64, Nothing}
    c_out         ::Union{Float64, Nothing}
    eff           ::Union{Float64, Nothing}
    eff_an        ::Union{Float64, Nothing}
    fluid         ::Union{Fluid, Nothing}
    membrane      ::Union{Membrane, Nothing}
    name          ::Union{String, Nothing}
    loss          ::Bool
    inv           ::Union{Float64, Nothing}
    p_out         ::Float64
    delta_p       ::Union{Float64, Nothing}
    U             ::Union{Float64, Nothing}
    pumping_power ::Union{Float64, Nothing}
    cost          ::Union{Float64, Nothing}
    pipe_flowrate ::Union{Float64, Nothing}
    flowrate      ::Union{Float64, Nothing}
    H             ::Union{Float64, Nothing}   # dimensionless (mass transport / surface)
    W             ::Union{Float64, Nothing}   # dimensionless (diffusion / surface)
    zeta          ::Union{Float64, Nothing}   # LM Sievert permeability ratio
    tau           ::Union{Float64, Nothing}   # dimensionless residence time
    alpha         ::Union{Float64, Nothing}   # MS permeability parameter [mol/m³]
    xi            ::Union{Float64, Nothing}   # MS saturation parameter (= alpha / c_in)
    J_perm        ::Union{Float64, Nothing}   # permeation flux [mol/m²/s]
    n_pipes       ::Union{Float64, Nothing}   # mirror of geometry.n_pipes
end

function Component(;
    geometry=nothing, c_in=nothing, c_out=nothing, eff=nothing,
    fluid=nothing, membrane=nothing, name=nothing, loss::Bool=false,
    inv=nothing, p_out::Float64=1e-15, delta_p=nothing, U=nothing,
    pumping_power=nothing, cost=nothing)
    n = geometry !== nothing ? geometry.n_pipes : nothing
    return Component(geometry, c_in, c_out, eff, nothing, fluid, membrane,
                     name, loss, inv, p_out, delta_p, U, pumping_power, cost,
                     nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing, n)
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
# These private functions encode the geometric and physical relationships that
# appear repeatedly in the permeation and inventory formulae.

"""
Cylindrical log-mean half-thickness of the membrane wall [m].

This is the correct characteristic length for 1-D radial diffusion through a
cylindrical shell:  denom = r_in · ln(r_out / r_in)

It appears in the denominator of every Fickian permeation-flux expression:
  J_diff = D_membrane / denom · K_S · (c_surface − c_equil)
"""
_wall_denom(d_hyd::Float64, dw::Float64)::Float64 =
    (d_hyd / 2) * log((d_hyd / 2 + dw) / (d_hyd / 2))
_wall_denom(comp::Component)::Float64 = _wall_denom(comp.fluid.d_Hyd, comp.membrane.dw)

"""Diffusion conductance of the membrane wall: G = D_m / wall_denom [m/s]."""
_diff_conductance(comp::Component)::Float64 = comp.membrane.D / _wall_denom(comp)

"""
LM Sievert permeability ratio ζ (zeta).

  ζ = D_m · K_S_membrane / (k_t · K_S_liquid · wall_denom)

Physical meaning: ratio of membrane diffusion resistance to fluid-side
mass-transport resistance.
  ζ >> 1  →  mass-transport limited
  ζ << 1  →  diffusion limited
"""
function _lm_zeta(comp::Component)::Float64
    return comp.membrane.D * comp.membrane.K_S /
           (comp.fluid.k_t * comp.fluid.Solubility * _wall_denom(comp))
end

"""
LM characteristic axial decay length [m⁻¹].

The local concentration decays along the pipe as c(L) ∝ exp(L_char · L),
where L_char < 0 for extraction.

  L_char = −ζ / (1 + ζ) · 4 k_t / (U₀ d_H)
"""
function _lm_L_char(comp::Component)::Float64
    z = _lm_zeta(comp)
    return -z / (1 + z) * 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
end

"""
MS permeability parameter α [mol/m³].

  α = (K_S · D_m / (4 k_t · wall_denom))² / K_H

Characterises the ratio of membrane permeance squared to the product of the
mass-transfer rate and the Henry constant. Used in the Lambert W analytical
solution for molten-salt extraction efficiency.
"""
function _ms_alpha(comp::Component)::Float64
    d = _wall_denom(comp)
    return (comp.membrane.K_S * comp.membrane.D / (4 * comp.fluid.k_t * d))^2 /
           comp.fluid.Solubility
end

"""
Safe evaluation of the principal Lambert W function W₀(exp(β)).

For very large β the argument exp(β) overflows Float64, so the asymptotic
approximation W(eˣ) ≈ x − ln(x) is used instead.
"""
function _lambertw_safe(beta_tau::Float64)::Float64
    beta_tau > log(floatmax(Float64)) && return beta_tau - log(beta_tau)
    return real(lambertw(complex(exp(beta_tau)), 0))
end

# ---------------------------------------------------------------------------
# update_attribute! override — propagates temperature changes
# ---------------------------------------------------------------------------

"""Override: a temperature update is forwarded to fluid and membrane."""
function update_attribute!(comp::Component, attr_name, new_value)::Nothing
    attr = attr_name isa Symbol ? attr_name : Symbol(attr_name)
    if attr === :T
        comp.fluid    !== nothing && (comp.fluid.T    = Float64(new_value))
        comp.membrane !== nothing && (comp.membrane.T = Float64(new_value))
        update_T_prop!(comp)
        return nothing
    end
    invoke(update_attribute!, Tuple{Any, Any, Any}, comp, attr_name, new_value)
    return nothing
end

# ---------------------------------------------------------------------------
# Utility methods
# ---------------------------------------------------------------------------

"""Re-evaluate temperature-dependent Arrhenius properties in fluid and membrane."""
function update_T_prop!(comp::Component)::Nothing
    comp.fluid    !== nothing && PipeSubclasses.update_T_prop!(comp.fluid)
    comp.membrane !== nothing && PipeSubclasses.update_T_prop!(comp.membrane)
    return nothing
end

"""
    friction_factor(comp, Re) -> f

Darcy-Weisbach friction factor. Uses the Blasius correlation (turbulent) or
the laminar result (64 / Re), with the transition at Re = 2300.
"""
friction_factor(::Component, Re::Float64)::Float64 =
    Re < 2300 ? 64.0 / Re : 0.316 / Re^0.25

"""
    get_pressure_drop!(comp) -> Δp [Pa]

Compute the Darcy-Weisbach pressure drop and store it in `comp.delta_p`.
"""
function get_pressure_drop!(comp::Component)::Float64
    Re_val = Correlations.Re(comp.fluid.rho, comp.fluid.U0, comp.geometry.D, comp.fluid.mu)
    f_val  = friction_factor(comp, Re_val)
    comp.delta_p = f_val * (comp.geometry.L / comp.geometry.D) *
                   (comp.fluid.rho * comp.fluid.U0^2) / 2.0
    return comp.delta_p
end

"""
    get_pipe_flowrate(comp) -> Q [m³/s]

Volumetric flow rate of a single pipe.
"""
function get_pipe_flowrate(comp::Component)::Float64
    comp.pipe_flowrate = comp.fluid.U0 * π * comp.fluid.d_Hyd^2 / 4.0
    return comp.pipe_flowrate
end

"""
    get_total_flowrate(comp) -> Q_total [m³/s]

Total volumetric flow rate summed over all `n_pipes` pipes.
"""
function get_total_flowrate(comp::Component)::Float64
    get_pipe_flowrate(comp)
    comp.flowrate = comp.pipe_flowrate * comp.geometry.n_pipes
    return comp.flowrate
end

"""
    get_pumping_power!(comp) -> P [W]

Compute pumping power for the component (all pipes in parallel).
Calls `get_pressure_drop!` if not already computed.
"""
function get_pumping_power!(comp::Component)::Float64
    comp.delta_p === nothing && get_pressure_drop!(comp)
    comp.pumping_power = comp.delta_p * get_pipe_flowrate(comp) * comp.geometry.n_pipes
    return comp.pumping_power
end

"""
    set_hydraulic_diameter!(comp) -> d_Hyd [m]

Derive the fluid hydraulic diameter from the component geometry's cross-section
(``D_h = 4 A / P``) and store it in `comp.fluid.d_Hyd`, overwriting any current
value. For a circular pipe this equals `comp.geometry.D`.
"""
function set_hydraulic_diameter!(comp::Component)::Float64
    comp.fluid.d_Hyd = hydraulic_diameter(comp.geometry)
    return comp.fluid.d_Hyd
end

"""
    define_component_volumes!(comp)

Copy fluid and wall volumes from `comp.geometry` into `comp.fluid.V` and
`comp.membrane.V`. If `comp.fluid.d_Hyd` is not already set, also derive it
from the geometry cross-section via [`set_hydraulic_diameter!`](@ref).
"""
function define_component_volumes!(comp::Component)::Nothing
    comp.fluid.V    = get_fluid_volume(comp.geometry)
    comp.membrane.V = get_solid_volume(comp.geometry)
    comp.fluid.d_Hyd === nothing && set_hydraulic_diameter!(comp)
    return nothing
end

"""
    estimate_cost!(comp; metal_cost=0.0, fluid_cost=0.0) -> cost [\$]

Estimate fabrication cost from volumetric cost rates [\$/m³].
Stores result in `comp.cost`.
"""
function estimate_cost!(comp::Component;
                        metal_cost::Float64=0.0, fluid_cost::Float64=0.0)::Float64
    comp.cost = (get_solid_volume(comp.geometry) * metal_cost +
                 get_fluid_volume(comp.geometry) * fluid_cost) * comp.geometry.n_pipes
    return comp.cost
end

"""
    T_leak(comp) -> mol/s

Tritium leakage rate from the component.
"""
T_leak(comp::Component)::Float64 = comp.c_in * comp.eff * get_pipe_flowrate(comp)

# ---------------------------------------------------------------------------
# Regime identification
# ---------------------------------------------------------------------------

"""
    get_regime(comp; print_var=false) -> String

Identify the dominant tritium transport regime: "Mass transport limited",
"Diffusion Limited", "Surface limited", or "Mixed regime".

Computes `k_t` if not already set.
"""
function get_regime(comp::Component; print_var::Bool=false)::String
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    if comp.fluid.MS
        return FusionCoolant.get_regime_ms(;
            k_d=comp.membrane.k_d, D=comp.membrane.D, dw=comp.membrane.dw,
            K_S=comp.membrane.K_S, c0=comp.c_in, k_t=comp.fluid.k_t,
            k_H=comp.fluid.Solubility, print_var=print_var)
    else
        return FusionCoolant.get_regime_lm(;
            D=comp.membrane.D, k_t=comp.fluid.k_t, K_S_S=comp.membrane.K_S,
            K_S_L=comp.fluid.Solubility, k_r=comp.membrane.k_r,
            dw=comp.membrane.dw, c0=comp.c_in, print_var=print_var)
    end
end

# ---------------------------------------------------------------------------
# Dimensionless parameters H and W
# ---------------------------------------------------------------------------

"""
    get_adimensionals!(comp)

Compute and store the dimensionless transport parameters `comp.H` and `comp.W`.

For molten salts:
  H = k_d / (k_t · K_H)    (surface-to-MT resistance ratio)
  W = 2 k_d δ √(c₀/K_H) / (D K_S)   (surface-to-diffusion ratio)

For liquid metals, the analogous `W_lm` and partition parameter are used.
"""
function get_adimensionals!(comp::Component)::Nothing
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    if comp.fluid.MS
        comp.H = FusionCoolant.H_ms(comp.fluid.k_t, comp.fluid.Solubility, comp.membrane.k_d)
        comp.W = FusionCoolant.W_ms(comp.membrane.k_d, comp.membrane.D, comp.membrane.dw,
                                   comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
    else
        comp.W = FusionCoolant.W_lm(comp.membrane.k_r, comp.membrane.D, comp.membrane.dw,
                                    comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
        comp.H = comp.W * FusionCoolant.partition_param_lm(comp.membrane.D, comp.fluid.k_t,
                                                           comp.membrane.K_S, comp.fluid.Solubility,
                                                           comp.membrane.dw)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Analytical efficiency
# ---------------------------------------------------------------------------

"""
    analytical_efficiency!(comp; p_out=1e-15)

Compute the closed-form extraction efficiency `comp.eff_an`.

**Molten salt** — uses the Lambert W solution derived in Humrickhouse et al.:
  - Mass-transport limited (ξ >> 1):  1 − exp(−τ)
  - Diffusion limited (√ξ << τ):     1 − (1 − τ√ξ)²
  - General case:                     Lambert W₀(exp(β − τ − 1))

**Liquid metal** — uses the exponential efficiency with the permeability
parameter ζ = membrane diffusion resistance / MT resistance:
  eff = (1 − exp(−τ ζ/(1+ζ))) · (1 − √(p_out/p_in))

Sets `comp.tau`, `comp.alpha`, `comp.xi` (MS) or `comp.zeta` (LM).
"""
function analytical_efficiency!(comp::Component; p_out::Float64=1e-15)::Nothing
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    comp.tau = 4.0 * comp.fluid.k_t * comp.geometry.L /
               (comp.fluid.U0 * comp.fluid.d_Hyd)

    if comp.fluid.MS
        comp.alpha = _ms_alpha(comp)
        comp.xi    = comp.alpha / comp.c_in
        p_in       = comp.c_in / comp.fluid.Solubility

        if comp.xi > 1e5
            # Mass-transport limited: membrane equilibrates fast, MT controls
            comp.eff_an = (1.0 - exp(-comp.tau)) * (1.0 - p_out / p_in)

        elseif comp.xi^0.5 < 1e-2 && comp.tau > 1.0 / comp.xi^0.5
            # Diffusion limited: linear concentration profile approximation
            comp.eff_an = (1.0 - (1.0 - comp.tau * comp.xi^0.5)^2) *
                          (1.0 - (p_out / p_in)^0.5)

        else
            # General case — Lambert W analytical solution
            e        = (comp.alpha * p_out * comp.fluid.Solubility)^0.5
            f        = e / comp.alpha
            delta    = (1.0 / comp.xi + 1.0 + 2 * f)^0.5
            beta     = delta + (1.0 + f) * log(delta - 1.0 - f)
            beta_tau = beta - comp.tau - 1.0

            if beta_tau > log(floatmax(Float64)) || p_out > 1e-5
                # Iterative fallback when exp(beta_tau) overflows or p_out is significant
                lo = min(p_out * comp.fluid.Solubility, comp.c_in)
                hi = max(p_out * comp.fluid.Solubility, comp.c_in)
                if abs(lo - hi) / hi < 1e-2
                    comp.eff_an = 1e-6
                    return nothing
                end
                function eq(cl)
                    lhs = (cl / comp.alpha + 1.0 + 2 * f)^0.5 +
                          (1.0 + f) * log(-f + ((cl / comp.alpha + 1.0 + 2 * f)^0.5 - 1.0))
                    return (lhs - (beta - comp.tau))^2
                end
                res = optimize(eq, lo, hi, Brent())
                comp.eff_an = 1.0 - Optim.minimizer(res) / comp.c_in
            else
                w = _lambertw_safe(beta_tau)
                comp.eff_an = 1.0 - comp.xi * (w^2 + 2w)
            end
        end

    else
        # Liquid metal — single exponential with permeability ratio ζ
        comp.zeta  = _lm_zeta(comp)
        p_in       = (comp.c_in / comp.fluid.Solubility)^2
        comp.eff_an = (1.0 - exp(-comp.tau * comp.zeta / (1.0 + comp.zeta))) *
                      (1.0 - (p_out / p_in)^0.5)
    end
    return nothing
end

"""
    use_analytical_efficiency!(comp; p_out=1e-15)

Convenience wrapper: compute `eff_an` and copy it into `comp.eff`.
"""
function use_analytical_efficiency!(comp::Component; p_out::Float64=1e-15)::Nothing
    analytical_efficiency!(comp; p_out=p_out)
    comp.eff = comp.eff_an
    return nothing
end

# ---------------------------------------------------------------------------
# Permeation flux helpers (internal)
# ---------------------------------------------------------------------------
# The flux helpers compute J_perm [mol/m²/s] for a given local bulk
# concentration c and set comp.J_perm. They return the wall concentration
# c_wl for use as the starting guess in the next axial step.
#
# Sign convention: J_perm < 0 means tritium leaves the fluid (extraction).
# NOTE: the MT+surface mixed branch (W < 0.1, 0.01 < H < 100) currently
# sets J_perm > 0 for the same physics, creating a sign inconsistency that
# causes the concentration to rise instead of fall in get_efficiency!.
# See trioma_jl_issues.md Issue A. The affected tests are marked @test_broken.

"""Pure diffusion-limited flux for molten salt."""
function _diffusion_flux_ms(comp::Component, c::Float64, p_out::Float64)::Float64
    return -_diff_conductance(comp) * comp.membrane.K_S *
           ((c / comp.fluid.Solubility)^0.5 - p_out^0.5)
end

"""Pure diffusion-limited flux for liquid metal."""
function _diffusion_flux_lm(comp::Component, c::Float64, p_out::Float64)::Float64
    return -_diff_conductance(comp) * comp.membrane.K_S *
           (c / comp.fluid.Solubility - p_out^0.5)
end

function _get_flux_ms!(comp::Component, c::Float64, W::Float64, H::Float64;
                       c_guess::Float64, p_out::Float64)::Float64
    G = _diff_conductance(comp)   # membrane diffusion conductance [m/s]

    if W > 10
        # Diffusion (W) dominates over surface kinetics
        if H / W > 1000
            # MT limited: membrane equilibrates fast, fluid-side MT controls
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H / W < 0.0001
            # Diffusion limited
            comp.J_perm = _diffusion_flux_ms(comp, c, p_out)
        else
            # Mixed MT + diffusion: find wall concentration by balancing J_mt = J_diff
            function eq_md(c_wl)
                J_mt   = 2 * comp.fluid.k_t * (c - c_wl)
                J_diff = G * comp.membrane.K_S * ((c_wl / comp.fluid.Solubility)^0.5 - p_out^0.5)
                return abs(J_diff - J_mt)
            end
            res = optimize(eq_md, 0.0, c, Brent())
            cwl = Optim.minimizer(res)
            comp.J_perm = -2 * comp.fluid.k_t * (c - cwl)
            return cwl
        end

    elseif W < 0.1
        # Surface kinetics (W) dominates over diffusion
        if H > 100
            # MT limited
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H < 0.01
            # Surface limited
            comp.J_perm = -comp.membrane.k_d * (c / comp.fluid.Solubility)
        else
            # Mixed MT + surface: balance J_mt = J_surf
            # NOTE: J_perm is set positive here — sign inconsistency vs other branches
            # (see trioma_jl_issues.md Issue A)
            function eq_ms(c_wl)
                J_mt   = 2 * comp.fluid.k_t * (c - c_wl)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                return abs(J_mt - J_surf)
            end
            res = optimize(eq_ms, 0.0, c, Brent())
            cwl = Optim.minimizer(res)
            comp.J_perm = 2 * comp.fluid.k_t * (c - cwl)
            return cwl
        end

    else
        # Intermediate W: diffusion–surface mixed regime
        if H / W > 1000
            # MT limited
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H / W < 0.0001
            # Balance diffusion vs surface at the outer wall
            function eq_sd(c_wl)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                J_diff = G * comp.membrane.K_S * ((c_wl / comp.fluid.Solubility)^0.5 - p_out^0.5)
                return abs(J_diff - J_surf)
            end
            res = optimize(eq_sd, 1e-14, c, Brent())
            cw  = Optim.minimizer(res)
            comp.J_perm = G * (comp.membrane.K_S * (cw / comp.fluid.Solubility)^0.5 - p_out^0.5)
            return cw
        else
            # Fully coupled: simultaneously balance MT, diffusion, and surface reactions
            function eq_full(cv)
                c_wl, c_ws = cv
                J_mt   = 2 * comp.fluid.k_t * (c - c_wl)
                J_d    = comp.membrane.k_d * (c_wl / comp.membrane.K_S) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_ws^2
                J_diff = G * (comp.membrane.K_S * c_ws - p_out^0.5)
                return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
            end
            res = optimize(eq_full, [0.0, 0.0], [c, c], [2c/3, c/3],
                           Fminbox(NelderMead()), Optim.Options(g_tol=1e-8))
            cwl = Optim.minimizer(res)[1]
            comp.J_perm = 2 * comp.fluid.k_t * (c - cwl)
            return cwl
        end
    end
    return c_guess   # unreachable: all branches above return explicitly
end

function _get_flux_lm!(comp::Component, c::Float64, W::Float64, H::Float64;
                       c_guess::Float64, p_out::Float64)::Float64
    G = _diff_conductance(comp)

    if W > 10
        if H / W > 1000
            comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
        elseif H / W < 0.0001
            comp.J_perm = _diffusion_flux_lm(comp, c, p_out)
        else
            # Mixed MT + diffusion
            function eq_md(c_wl)
                J_mt   = comp.fluid.k_t * (c - c_wl)
                J_diff = G * comp.membrane.K_S * (c_wl / comp.fluid.Solubility - p_out^0.5)
                return abs(J_diff - J_mt)
            end
            res = optimize(eq_md, 0.0, c, Brent())
            cwl = Optim.minimizer(res)
            comp.J_perm = -comp.fluid.k_t * (c - cwl)
            return cwl
        end

    elseif W < 0.1
        if H > 100
            comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
        elseif H < 0.01
            comp.J_perm = -comp.membrane.k_d * (c / comp.fluid.Solubility)
        else
            # Mixed MT + surface
            function eq_ms(c_wl)
                J_mt   = comp.fluid.k_t * (c - c_wl - p_out^0.5 * comp.fluid.Solubility)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                return abs(J_mt - J_surf)
            end
            res = optimize(eq_ms, 0.0, c, Brent())
            cwl = Optim.minimizer(res)
            comp.J_perm = comp.fluid.k_t * (c - cwl - p_out * comp.fluid.Solubility)
            return cwl
        end

    else
        if H / W < 0.0001
            comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
        elseif H / W > 1000
            function eq_sd(c_wl)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                J_diff = G * comp.membrane.K_S * (c_wl / comp.fluid.Solubility - p_out^0.5)
                return abs(J_diff - J_surf)
            end
            res = optimize(eq_sd, 1e-14, c, Brent())
            cw  = Optim.minimizer(res)
            comp.J_perm = G * comp.membrane.K_S * (cw / comp.fluid.Solubility - p_out^0.5)
            return cw
        else
            # Fully coupled MT + diffusion + surface
            function eq_full(cv)
                c_wl, c_ws = cv
                J_mt   = comp.fluid.k_t * (c - c_wl)
                J_d    = comp.membrane.k_d * (c_wl / comp.membrane.K_S) -
                         comp.membrane.k_d * comp.membrane.K_S * c_ws^2
                J_diff = G * (comp.membrane.K_S * c_ws - p_out^0.5)
                return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
            end
            ub  = c * (1 + 1e-4)
            res = optimize(eq_full, [0.0, 0.0], [ub, ub], [2c/3, c/3],
                           Fminbox(NelderMead()), Optim.Options(g_tol=1e-8))
            cwl = Optim.minimizer(res)[1]
            comp.J_perm = comp.fluid.k_t * (c - cwl)
            return cwl
        end
    end
    return c_guess   # unreachable: all branches above return explicitly
end

# ---------------------------------------------------------------------------
# Public flux interface
# ---------------------------------------------------------------------------

"""
    get_flux!(comp, c; c_guess=1e-9, p_out=1e-15) -> J_perm

Compute the local permeation flux `J_perm` [mol/m²/s] at bulk concentration `c`,
store it in `comp.J_perm`, and return it.

Dispatches on `comp.fluid.MS` and the dimensionless regime parameters H and W.
Calls `get_adimensionals!` internally. Sign convention: J_perm < 0 means tritium
leaves the fluid (extraction toward vacuum side).
"""
function get_flux!(comp::Component, c::Float64;
                   c_guess::Float64=1e-9, p_out::Float64=1e-15)::Float64
    get_adimensionals!(comp)
    W = comp.W; H = comp.H
    if comp.fluid.MS
        _get_flux_ms!(comp, c, W, H; c_guess=c_guess, p_out=p_out)
    else
        _get_flux_lm!(comp, c, W, H; c_guess=c_guess, p_out=p_out)
    end
    return comp.J_perm
end

# ---------------------------------------------------------------------------
# Numerical efficiency (axial integration)
# ---------------------------------------------------------------------------

"""
    get_efficiency!(comp; c_guess=nothing, p_out=1e-15)

Numerically integrate the tritium concentration profile along the pipe length
and store the extraction efficiency in `comp.eff`.

Uses 100 uniform axial steps. At each step the local flux is computed and the
concentration is advanced using the plug-flow mass balance:

  Δc = f_H₂ · J_perm · (4 / U₀ / d_H) · ΔL

where `f_H₂ = 0.5` for molten-salt (tritium dissolves as T₂, so each permeating
T₂ molecule removes 2 dissolved T atoms) and `f_H₂ = 1.0` for liquid-metal
(atomic tritium).
"""
function get_efficiency!(comp::Component;
                         c_guess::Union{Float64,Nothing}=nothing,
                         p_out::Float64=1e-15)::Nothing
    if comp.c_in == 0
        comp.c_out = 0.0; comp.eff = 0.0; return nothing
    end
    L_vec = range(0, comp.geometry.L, length=100)
    dl    = step(L_vec)
    c_vec = zeros(Float64, length(L_vec))
    # Stoichiometric factor: T₂ permeation removes 2 T atoms from an MS fluid
    f_H2  = comp.fluid.MS ? 0.5 : 1.0
    cg    = c_guess !== nothing ? c_guess : Float64(comp.c_in)

    # Compute H and W once outside the loop; use private helpers directly so
    # that the returned c_wl (inner wall concentration) serves as the warm-start
    # guess for the next axial step without interfering with get_flux!'s return value.
    get_adimensionals!(comp)
    W = comp.W; H = comp.H

    for i in eachindex(L_vec)
        if i == 1
            c_vec[1] = Float64(comp.c_in)
        else
            c_vec[i] = c_vec[i-1] + f_H2 * comp.J_perm * 4 * dl /
                       (comp.fluid.U0 * comp.fluid.d_Hyd)
        end
        cg = if comp.fluid.MS
            _get_flux_ms!(comp, c_vec[i], W, H; c_guess=cg, p_out=p_out)
        else
            _get_flux_lm!(comp, c_vec[i], W, H; c_guess=cg, p_out=p_out)
        end
    end
    comp.eff = (comp.c_in - c_vec[end]) / comp.c_in
    return nothing
end

# ---------------------------------------------------------------------------
# Outlet concentration with recirculation
# ---------------------------------------------------------------------------

"""
    outlet_c_comp!(comp) -> c_out [mol/m³]

Compute the outlet concentration from the current extraction efficiency,
accounting for recirculation.

- `recirculation == 0`:  c_out = c_in · (1 − eff)  (simple pass-through)
- `recirculation > 0`:   fraction `r` of the outlet is recycled back to the
  inlet; iterates to self-consistency
- `recirculation < 0`:   fraction `|r|` is bypassed; c_out is a blend of
  extracted and bypassed streams
"""
function outlet_c_comp!(comp::Component)::Float64
    r = comp.fluid.recirculation
    if r == 0.0
        comp.c_out = comp.c_in * (1.0 - comp.eff)

    elseif r > 0.0
        comp.c_in == 0.0 && error("Inlet concentration is zero with positive recirculation")
        c0 = comp.c_in
        c  = c0
        tol = 1e-6; err = 1.0
        while err > tol
            c1 = c
            update_attribute!(comp, "c_in", c)
            analytical_efficiency!(comp)
            update_attribute!(comp, "eff", comp.eff_an)
            comp.c_out = comp.c_in * (1.0 - comp.eff)
            c = (comp.c_out * r + c0) / (r + 1.0)
            err = abs((c - c1) / c)
        end

    elseif r < 0.0
        r <= -1.0 && error("Bypass (negative recirculation) ≤ -1 is invalid")
        comp.c_in == 0.0 && error("Inlet concentration is zero with bypass recirculation")
        comp.c_out = comp.c_in * (1.0 - comp.eff) * (1.0 + r) + comp.c_in * (-r)
    end
    return comp.c_out
end

# ---------------------------------------------------------------------------
# Global heat-exchanger coefficient
# ---------------------------------------------------------------------------

"""
    get_global_HX_coeff!(comp; R_conv_sec=0.0)

Compute the overall heat-transfer coefficient U [W/m²/K] for the tube-wall.

Combines the convective (fluid-side) and conductive (wall) resistances:
  1/U = 1/h_prim + R_cond + R_conv_sec

`R_conv_sec` is the optional secondary-side convective resistance [m²K/W].
Stores h_prim in `comp.fluid.h_coeff` and U in `comp.U`.
"""
function get_global_HX_coeff!(comp::Component; R_conv_sec::Float64=0.0)::Nothing
    R_cond = log((comp.fluid.d_Hyd + comp.membrane.dw) / comp.fluid.d_Hyd) /
             (2π * comp.membrane.k)
    Re_val = Correlations.Re(comp.fluid.rho, comp.fluid.U0, comp.fluid.d_Hyd, comp.fluid.mu)
    Pr_val = Correlations.Pr(comp.fluid.cp, comp.fluid.mu, comp.fluid.k)

    if comp.geometry.turbulator === nothing
        h_prim = Correlations.get_h_from_Nu(
            Correlations.Nu_DittusBoelter(Re_val, Pr_val), comp.fluid.k, comp.fluid.d_Hyd)
    elseif isa(comp.geometry.turbulator, WireCoil)
        h_prim = PipeSubclasses.h_t_correlation(comp.geometry.turbulator;
                     Re=Re_val, Pr=Pr_val, d_hyd=comp.fluid.d_Hyd, k=comp.fluid.k)
    elseif isa(comp.geometry.turbulator, CustomTurbulator)
        h_prim = PipeSubclasses.h_t_correlation(comp.geometry.turbulator;
                     Re=Re_val, Pr=Pr_val, d_hyd=comp.fluid.d_Hyd, k=comp.fluid.k)
    else
        error("Unknown turbulator type: $(comp.geometry.turbulator.turbulator_type)")
    end
    comp.fluid.h_coeff = h_prim
    comp.U = 1.0 / (1.0 / h_prim + R_cond + R_conv_sec)
    return nothing
end

# ---------------------------------------------------------------------------
# Inventory helpers (internal)
# ---------------------------------------------------------------------------
# Both analytical and numerical approaches compute a 2-D integral over the
# pipe length L and the radial coordinate r through the membrane wall.
# The concentration profile c_m(r,L) is separable:
#
#   LM: c_m(r,L) = (c_w(L) − c_ext) · (−log(r/r_out)/log(r_out/r_in)) + c_ext
#   MS: same radial profile; c_wl(L) from the Lambert W solution
#
# where r_out = r_in + dw = D/2 + dw.

# Closed-form radial integral of the log-profile shape function:
# ∫_r_in^r_out (−log(r/r_out) / log(r_out/r_in)) · 2πr dr
# Antiderivative: ifun(r) = r²/4 · (2·log(r/r_out) − 1)  →  d/dr[ifun] = r·log(r/r_out)
_radial_log_integral(r_in::Float64, r_out::Float64)::Float64 =
    let ifun(r) = r^2/4 * (2*log(r/r_out) - 1)
        2π / log(r_out/r_in) * (ifun(r_in) - ifun(r_out))
    end

# ---------------------------------------------------------------------------
# Solid inventory (membrane)
# ---------------------------------------------------------------------------

"""
    get_solid_inventory!(comp; p_out=0.0, flag_an=false) -> inv [mol]

Membrane (solid) tritium inventory. Default: 1-D numerical integration (single
`quadgk` over the pipe length after analytically evaluating the radial integral).
Pass `flag_an=true` to use the closed-form formula instead.

The concentration profile through the membrane wall is:
    c_m(r,L) = (−log(r/r_out)/log(r_out/r_in)) · (c_wl(L) − c_ext) + c_ext
The r-integral of c_m(r,L)·2πr has the closed form R_shape·(c_wl−c_ext) + c_ext·π(r_out²−r_in²),
reducing the 2-D problem to a single 1-D quadgk over L.

Stores the result in `comp.membrane.inv`.
"""
function get_solid_inventory!(comp::Component;
                               p_out::Float64=0.0, flag_an::Bool=false)::Float64
    flag_an && return analytical_solid_inventory!(comp; p_out=p_out)
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    r_in  = comp.fluid.d_Hyd / 2
    r_out = r_in + comp.membrane.dw
    c_ext = p_out^0.5 * comp.membrane.K_S

    # Closed-form radial integral — constants hoisted outside quadgk
    R_shape    = _radial_log_integral(r_in, r_out)
    r_cross    = π * (r_out^2 - r_in^2)   # cross-sectional annular area

    # Wall-concentration closure: hoist all L-independent constants
    cwl = if comp.fluid.MS
        α         = _ms_alpha(comp)
        xi_l      = α / comp.c_in
        beta      = (1/xi_l + 1)^0.5 + log((1/xi_l + 1)^0.5 - 1)
        tau_coef  = 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
        wl_scale  = (α / comp.fluid.Solubility)^0.5 * comp.membrane.K_S
        L -> wl_scale * _lambertw_safe(beta - tau_coef * L - 1.0) + c_ext
    else
        z         = _lm_zeta(comp)
        L_ch      = _lm_L_char(comp)
        A_lm      = comp.c_in * comp.membrane.K_S / (comp.fluid.Solubility * (z + 1))
        L -> A_lm * exp(L_ch * L) + c_ext
    end

    result, _ = quadgk(L -> R_shape * (cwl(L) - c_ext) + c_ext * r_cross,
                       0.0, comp.geometry.L)
    comp.membrane.inv = result * comp.geometry.n_pipes
    isnan(comp.membrane.inv) && println("Error: Inventory calculation failed")
    return comp.membrane.inv
end

"""
    analytical_solid_inventory!(comp; p_out=0.0) -> inv [mol]

Closed-form membrane inventory (LM only). For molten salt this function
delegates to the numerical integrator because the MS analytical formula has
a known sign error (see trioma_jl_issues.md Bug 1).

For LM, integrates the exact exponential concentration profile analytically
using the antiderivative:  ifun(r) = r²/4 · (2·ln(r/r_out) − 1)

The full inventory integral is:
  inv = K · [ifun(r_out) − ifun(r_in)] + L_pipe · c_ext · π · (r_out² − r_in²)

where K is the axial prefactor that folds in the exponential decay along L.
"""
function analytical_solid_inventory!(comp::Component; p_out::Float64=0.0)::Float64
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    if !comp.fluid.MS
        r_in  = comp.geometry.D / 2
        r_out = r_in + comp.geometry.dw
        z     = _lm_zeta(comp)
        L_ch  = _lm_L_char(comp)
        # K is the radial × azimuthal prefactor for the varying part of c_m
        K = -2π * comp.c_in * comp.membrane.K_S /
            (comp.fluid.Solubility * (z + 1) * log(r_out / r_in))
        # Fold in the axial integral ∫₀^L exp(L_ch·l) dl
        K *= (exp(L_ch * comp.geometry.L) - 1) / L_ch
        ifun(r) = r^2 / 4 * (2 * log(r / r_out) - 1)
        c_ext   = p_out^0.5 * comp.membrane.K_S
        integral = K * (ifun(r_out) - ifun(r_in)) +
                   comp.geometry.L * c_ext * π * (r_out^2 - r_in^2)
        comp.membrane.inv = integral
        return integral

    else
        # MS analytical formula has a sign error in the ifun expression; use numerical
        # See trioma_jl_issues.md Bug 1 for the full derivation of the correct formula
        return get_solid_inventory!(comp; flag_an=false, p_out=p_out)
    end
end

# ---------------------------------------------------------------------------
# Fluid inventory
# ---------------------------------------------------------------------------

"""
    get_fluid_inventory!(comp; flag_an=false, p_out=0.0) -> inv [mol]

Tritium inventory in the fluid phase. Default: numerical integration along L.
Pass `flag_an=true` to use the analytical formula (LM only; MS falls back to
numerical).

The fluid concentration profile c(L) follows from the plug-flow balance, so
the inventory is: π r_in² · ∫₀^L c(l) dl.
"""
function get_fluid_inventory!(comp::Component;
                               flag_an::Bool=false, p_out::Float64=0.0)::Float64
    flag_an && return analytical_fluid_inventory!(comp; p_out=p_out)
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    r_in = comp.fluid.d_Hyd / 2

    function integrand(L)
        if !comp.fluid.MS
            L_ch  = _lm_L_char(comp)
            c_ext = p_out^0.5 * comp.fluid.Solubility
            return (comp.c_in - c_ext) * exp(L_ch * L) + c_ext
        else
            α    = _ms_alpha(comp)
            xi   = α / comp.c_in
            tau  = 4 * comp.fluid.k_t * L / (comp.fluid.U0 * comp.fluid.d_Hyd)
            beta = (1/xi + 1)^0.5 + log((1/xi + 1)^0.5 - 1)
            w    = _lambertw_safe(beta - tau - 1.0)
            return α * (w^2 + 2w)
        end
    end

    result, _ = quadgk(integrand, 0.0, comp.geometry.L)
    comp.fluid.inv = result * π * r_in^2 * comp.geometry.n_pipes
    return comp.fluid.inv
end

"""
    analytical_fluid_inventory!(comp; p_out=0.0) -> inv [mol]

Closed-form fluid inventory (LM only). For molten salt prints a notice and
delegates to the numerical integrator.

For LM the concentration profile is c(L) = (c_in − c_ext)·exp(L_char·L) + c_ext,
and its integral ∫₀^L c(l) dl can be computed analytically.
"""
function analytical_fluid_inventory!(comp::Component; p_out::Float64=0.0)::Float64
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    if !comp.fluid.MS
        L_ch  = _lm_L_char(comp)
        c_ext = p_out^0.5 * comp.fluid.Solubility
        K     = (comp.c_in - c_ext) * (exp(L_ch * comp.geometry.L) - 1) / L_ch
        r_in  = comp.fluid.d_Hyd / 2
        integral = K * π * r_in^2 + c_ext * π * r_in^2 * comp.geometry.L
        comp.fluid.inv = integral * comp.geometry.n_pipes
        return integral
    else
        println("MS fluid integration is done numerically")
        return get_fluid_inventory!(comp; flag_an=false, p_out=p_out)
    end
end

"""
    get_inventory!(comp; flag_an=true, p_out=0.0)

Compute both solid and fluid inventories and store their sum in `comp.inv` [mol].
"""
function get_inventory!(comp::Component; flag_an::Bool=true, p_out::Float64=0.0)::Nothing
    get_solid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    get_fluid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    comp.inv = comp.fluid.inv + comp.membrane.inv
    return nothing
end

end # module PAVPipe
