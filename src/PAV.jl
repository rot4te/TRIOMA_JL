"""
PAV.jl

Component — the core permeation-against-vacuum (PAV) / heat-exchanger type
in TRIOMA. Contains analytical and numerical methods for:
  - transport-regime identification
  - analytical and numerical tritium extraction efficiency
  - permeation flux (all regimes)
  - tritium inventory in fluid and membrane
  - heat-exchanger splitting for temperature discretization
"""

module PAVModule

using ..TriomaTypes
import ..TriomaTypes: update_attribute!
using ..PipeSubclasses
using ..Correlations
using ..MoltenSalts
using ..LiquidMetals
using LambertW
using QuadGK
using Optim

export Component,
       get_regime, get_pipe_flowrate, get_total_flowrate,
       define_component_volumes!, get_adimensionals!,
       use_analytical_efficiency!, analytical_efficiency!,
       get_efficiency!, get_flux!, outlet_c_comp!,
       get_global_HX_coeff!, friction_factor, get_pressure_drop!,
       get_pumping_power!, estimate_cost!,
       get_solid_inventory!, get_fluid_inventory!, get_inventory!,
       analytical_solid_inventory!, analytical_fluid_inventory!,
       T_leak, update_T_prop!

# ---------------------------------------------------------------------------
# Component
# ---------------------------------------------------------------------------

"""
A single permeation or heat-exchanger pipe component.

All mutable fields; keyword constructor provided.
"""
mutable struct Component <: TriomaClass
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
    H             ::Union{Float64, Nothing}   # dimensionless (mass transport/surface)
    W             ::Union{Float64, Nothing}   # dimensionless (diffusion/surface)
    zeta          ::Union{Float64, Nothing}   # LM parameter
    tau           ::Union{Float64, Nothing}   # dimensionless time
    alpha         ::Union{Float64, Nothing}   # MS surface parameter
    xi            ::Union{Float64, Nothing}   # MS extraction parameter
    J_perm        ::Union{Float64, Nothing}   # permeation flux [mol/m²/s]
    n_pipes       ::Union{Float64, Nothing}   # (mirror of geometry.n_pipes)
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
# update_attribute! override — handles T propagation
# ---------------------------------------------------------------------------

function update_attribute!(comp::Component, attr_name, new_value)
    attr = attr_name isa Symbol ? attr_name : Symbol(attr_name)
    if attr === :T
        comp.fluid    !== nothing && (comp.fluid.T    = Float64(new_value))
        comp.membrane !== nothing && (comp.membrane.T = Float64(new_value))
        update_T_prop!(comp)
        return
    end
    invoke(update_attribute!, Tuple{TriomaClass, Any, Any}, comp, attr_name, new_value)
end

# ---------------------------------------------------------------------------
# Auxiliary methods
# ---------------------------------------------------------------------------

"""Update temperature-dependent properties in fluid and membrane."""
function update_T_prop!(comp::Component)
    comp.fluid    !== nothing && PipeSubclasses.update_T_prop!(comp.fluid)
    comp.membrane !== nothing && PipeSubclasses.update_T_prop!(comp.membrane)
end

"""Laminar/turbulent Darcy-Weisbach friction factor."""
function friction_factor(comp::Component, Re::Float64)::Float64
    Re < 2300 ? 64.0 / Re : 0.316 / Re^0.25
end

"""Pressure drop across the component [Pa]."""
function get_pressure_drop!(comp::Component)::Float64
    Re_val = Correlations.Re(comp.fluid.rho, comp.fluid.U0, comp.geometry.D, comp.fluid.mu)
    f_val  = friction_factor(comp, Re_val)
    comp.delta_p = f_val * (comp.geometry.L / comp.geometry.D) *
                   (comp.fluid.rho * comp.fluid.U0^2) / 2.0
    return comp.delta_p
end

"""Volumetric flow rate of one pipe [m³/s]."""
function get_pipe_flowrate(comp::Component)::Float64
    comp.pipe_flowrate = comp.fluid.U0 * π * comp.fluid.d_Hyd^2 / 4.0
    return comp.pipe_flowrate
end

"""Total volumetric flow rate [m³/s] summed over all pipes."""
function get_total_flowrate(comp::Component)::Float64
    get_pipe_flowrate(comp)
    comp.flowrate = comp.pipe_flowrate * comp.geometry.n_pipes
    return comp.flowrate
end

"""Pumping power [W] for the component."""
function get_pumping_power!(comp::Component)::Float64
    comp.delta_p === nothing && get_pressure_drop!(comp)
    comp.pumping_power = comp.delta_p * get_pipe_flowrate(comp) * comp.geometry.n_pipes
    return comp.pumping_power
end

"""Assign fluid and membrane volumes from geometry."""
function define_component_volumes!(comp::Component)
    comp.fluid.V    = get_fluid_volume(comp.geometry)
    comp.membrane.V = get_solid_volume(comp.geometry)
end

"""Estimated fabrication cost [\$]."""
function estimate_cost!(comp::Component;
                        metal_cost::Float64=0.0, fluid_cost::Float64=0.0)::Float64
    V_solid = get_solid_volume(comp.geometry)
    V_fluid = get_fluid_volume(comp.geometry)
    comp.cost = (V_solid * metal_cost + V_fluid * fluid_cost) * comp.geometry.n_pipes
    return comp.cost
end

"""Tritium leakage [mol/s]."""
function T_leak(comp::Component)::Float64
    return comp.c_in * comp.eff * get_pipe_flowrate(comp)
end

# ---------------------------------------------------------------------------
# Transport-regime identification
# ---------------------------------------------------------------------------

"""
    get_regime(comp) -> String

Identify the dominant tritium transport regime (mass transport / diffusion /
surface / mixed) for the component.
"""
function get_regime(comp::Component; print_var::Bool=false)::String
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    if comp.fluid.MS
        return MoltenSalts.get_regime_ms(;
            k_d=comp.membrane.k_d, D=comp.membrane.D, thick=comp.membrane.thick,
            K_S=comp.membrane.K_S, c0=comp.c_in, k_t=comp.fluid.k_t,
            k_H=comp.fluid.Solubility, print_var=print_var)
    else
        return LiquidMetals.get_regime_lm(;
            D=comp.membrane.D, k_t=comp.fluid.k_t, K_S_S=comp.membrane.K_S,
            K_S_L=comp.fluid.Solubility, k_r=comp.membrane.k_r,
            thick=comp.membrane.thick, c0=comp.c_in, print_var=print_var)
    end
end

# ---------------------------------------------------------------------------
# Dimensionless parameters H and W
# ---------------------------------------------------------------------------

"""
    get_adimensionals!(comp)

Compute and store the dimensionless transport parameters H and W.
"""
function get_adimensionals!(comp::Component)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    if comp.fluid.MS
        comp.H = MoltenSalts.H_ms(comp.fluid.k_t, comp.fluid.Solubility, comp.membrane.k_d)
        comp.W = MoltenSalts.W_ms(comp.membrane.k_d, comp.membrane.D, comp.membrane.thick,
                                   comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
    else
        comp.W = LiquidMetals.W_lm(comp.membrane.k_r, comp.membrane.D, comp.membrane.thick,
                                    comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
        comp.H = comp.W * LiquidMetals.partition_param_lm(comp.membrane.D, comp.fluid.k_t,
                                                           comp.membrane.K_S, comp.fluid.Solubility,
                                                           comp.membrane.thick)
    end
end

# ---------------------------------------------------------------------------
# Analytical efficiency
# ---------------------------------------------------------------------------

"""
    analytical_efficiency!(comp; p_out=1e-15)

Compute the analytical extraction efficiency `eff_an` and store it in `comp`.
Handles molten-salt (Lambert W) and liquid-metal regimes.
"""
function analytical_efficiency!(comp::Component; p_out::Float64=1e-15)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    comp.tau = 4.0 * comp.fluid.k_t * comp.geometry.L /
               (comp.fluid.U0 * comp.fluid.d_Hyd)

    if comp.fluid.MS
        # --- Molten salt ---
        comp.alpha = (1.0 / comp.fluid.Solubility *
            (0.5 * comp.membrane.K_S * comp.membrane.D /
             (comp.fluid.k_t * comp.fluid.d_Hyd *
              log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))))^2
        comp.xi  = comp.alpha / comp.c_in
        p_in     = comp.c_in / comp.fluid.Solubility

        if comp.xi > 1e5
            # Mass-transport limited
            corr_p = 1.0 - (p_out / p_in)
            comp.eff_an = (1.0 - exp(-comp.tau)) * corr_p

        elseif comp.xi^0.5 < 1e-2 && comp.tau > 1.0 / comp.xi^0.5
            # Diffusion-limited linear approximation
            corr_p = 1.0 - (p_out / p_in)^0.5
            comp.eff_an = (1.0 - (1.0 - comp.tau * comp.xi^0.5)^2) * corr_p

        else
            # General case — Lambert W
            e    = (comp.alpha * p_out * comp.fluid.Solubility)^0.5
            f    = e / comp.alpha
            delta = (1.0 / comp.xi + 1.0 + 2 * f)^0.5
            beta  = delta + (1.0 + f) * log(delta - 1.0 - f)
            max_exp   = log(floatmax(Float64))
            beta_tau  = beta - comp.tau - 1.0

            if beta_tau > max_exp || p_out > 1e-5
                # Solve iteratively (no Lambert W approximation available here)
                function eq(cl)
                    left  = (cl / comp.alpha + 1.0 + 2 * f)^0.5 +
                            (1.0 + f) * log(-f + ((cl / comp.alpha + 1.0 + 2 * f)^0.5 - 1.0))
                    right = beta - comp.tau
                    return (left - right)^2
                end
                lo1 = min(p_out * comp.fluid.Solubility, comp.c_in)
                hi1 = max(p_out * comp.fluid.Solubility, comp.c_in)
                if abs(p_out * comp.fluid.Solubility - comp.c_in) / comp.c_in < 1e-2
                    comp.eff_an = 1e-6
                    return
                end
                res = optimize(eq, lo1, hi1, Brent())
                comp.eff_an = 1.0 - (Optim.minimizer(res) / comp.c_in)
            else
                z = exp(beta_tau)
                w = lambertw(complex(z), 0)
                comp.eff_an = 1.0 - comp.xi * (real(w)^2 + 2 * real(w))
            end
        end

    else
        # --- Liquid metal ---
        comp.zeta = (2.0 * comp.membrane.K_S * comp.membrane.D) /
                    (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
                     log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
        p_in   = (comp.c_in / comp.fluid.Solubility)^2
        corr_p = 1.0 - (p_out / p_in)^0.5
        comp.eff_an = (1.0 - exp(-comp.tau * comp.zeta / (1.0 + comp.zeta))) * corr_p
    end
end

"""
    use_analytical_efficiency!(comp; p_out=1e-15)

Compute `eff_an` and copy it into `eff`.
"""
function use_analytical_efficiency!(comp::Component; p_out::Float64=1e-15)
    analytical_efficiency!(comp; p_out=p_out)
    comp.eff = comp.eff_an
end

# ---------------------------------------------------------------------------
# Numerical flux and efficiency
# ---------------------------------------------------------------------------

"""
    get_flux!(comp, c; c_guess=1e-9, p_out=1e-15) -> Float64

Compute the permeation flux J_perm [mol/m²/s] at local concentration `c`.
Returns the wall concentration for use as the next iteration's guess.
Dispatches on fluid.MS and transport regime (H, W).
"""
function get_flux!(comp::Component, c::Float64;
                   c_guess::Float64=1e-9, p_out::Float64=1e-15)::Float64
    get_adimensionals!(comp)
    W = comp.W; H = comp.H

    if comp.fluid.MS
        return _get_flux_ms!(comp, c, W, H; c_guess=c_guess, p_out=p_out)
    else
        return _get_flux_lm!(comp, c, W, H; c_guess=c_guess, p_out=p_out)
    end
end

# ---- Internal MS flux helpers ----
function _diffusion_flux_ms(comp, c, p_out)
    log_term = log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))
    denom    = comp.fluid.d_Hyd / 2 * log_term
    return -(comp.membrane.D / denom * comp.membrane.K_S *
             ((c / comp.fluid.Solubility)^0.5 - p_out^0.5))
end

function _get_flux_ms!(comp, c, W, H; c_guess, p_out)
    if W > 10
        if H / W > 1000
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H / W < 0.0001
            comp.J_perm = _diffusion_flux_ms(comp, c, p_out)
        else
            # Mixed: mass transport + diffusion
            log_term = log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))
            denom    = comp.fluid.d_Hyd / 2 * log_term
            function eq_md(c_wl)
                J_mt   = 2 * comp.fluid.k_t * (c - c_wl)
                J_diff = comp.membrane.D / denom * comp.membrane.K_S *
                         ((c_wl / comp.fluid.Solubility)^0.5 - p_out^0.5)
                return abs(J_diff - J_mt)
            end
            res = optimize(eq_md, 0.0, c, Brent())
            cwl = Optim.minimizer(res)
            comp.J_perm = -2 * comp.fluid.k_t * (c - cwl)
            return cwl
        end
    elseif W < 0.1
        if H > 100
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H < 0.01
            comp.J_perm = -comp.membrane.k_d * (c / comp.fluid.Solubility)
        else
            # Mixed: mass transport + surface
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
        # Mixed diffusion-surface
        if H / W > 1000
            comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
        elseif H / W < 0.0001
            log_term = log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))
            denom    = comp.fluid.d_Hyd / 2 * log_term
            function eq_sd(c_wl)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                J_diff = comp.membrane.D / denom *
                         comp.membrane.K_S * ((c_wl / comp.fluid.Solubility)^0.5 - p_out^0.5)
                return abs(J_diff - J_surf)
            end
            res = optimize(eq_sd, 1e-14, c, Brent())
            cw = Optim.minimizer(res)
            log_term2 = log((comp.fluid.d_Hyd/2 + comp.membrane.thick)/(comp.fluid.d_Hyd/2))
            comp.J_perm = comp.membrane.D / (comp.fluid.d_Hyd/2 * log_term2) *
                          (comp.membrane.K_S * (cw / comp.fluid.Solubility)^0.5 - p_out^0.5)
            return cw
        else
            # Fully coupled: mass transport + diffusion + surface
            function eq_full(cv)
                c_wl, c_ws = cv
                J_mt   = 2 * comp.fluid.k_t * (c - c_wl)
                log_t  = log((comp.fluid.d_Hyd/2 + comp.membrane.thick)/(comp.fluid.d_Hyd/2))
                J_d    = comp.membrane.k_d * (c_wl / comp.membrane.K_S) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_ws^2
                J_diff = comp.membrane.D / (comp.fluid.d_Hyd/2 * log_t) *
                         (comp.membrane.K_S * c_ws - p_out^0.5)
                return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
            end
            res = optimize(eq_full, [0.0, 0.0], [c, c], [2c/3, c/3], Fminbox(NelderMead()),
                           Optim.Options(g_tol=1e-8))
            cwl = Optim.minimizer(res)[1]
            comp.J_perm = 2 * comp.fluid.k_t * (c - cwl)
            return cwl
        end
    end
    return c_guess   # fallback
end

# ---- Internal LM flux helpers ----
function _diffusion_flux_lm(comp, c, p_out)
    log_term = log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))
    denom    = comp.fluid.d_Hyd / 2 * log_term
    return -(comp.membrane.D / denom * (comp.membrane.K_S *
             (c / comp.fluid.Solubility - p_out^0.5)))
end

function _get_flux_lm!(comp, c, W, H; c_guess, p_out)
    if W > 10
        if H / W > 1000
            comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
        elseif H / W < 0.0001
            comp.J_perm = _diffusion_flux_lm(comp, c, p_out)
        else
            log_term = log((comp.fluid.d_Hyd/2 + comp.membrane.thick)/(comp.fluid.d_Hyd/2))
            denom    = comp.fluid.d_Hyd / 2 * log_term
            function eq_md(c_wl)
                J_mt   = comp.fluid.k_t * (c - c_wl)
                J_diff = comp.membrane.D / denom *
                         (comp.membrane.K_S * (c_wl / comp.fluid.Solubility - p_out^0.5))
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
            log_term = log((comp.fluid.d_Hyd/2 + comp.membrane.thick)/(comp.fluid.d_Hyd/2))
            denom    = comp.fluid.d_Hyd / 2 * log_term
            function eq_sd(c_wl)
                J_surf = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                         comp.membrane.k_d * comp.membrane.K_S^2 * c_wl^2
                J_diff = comp.membrane.D / denom *
                         (comp.membrane.K_S * (c_wl / comp.fluid.Solubility - p_out^0.5))
                return abs(J_diff - J_surf)
            end
            res = optimize(eq_sd, 1e-14, c, Brent())
            cw = Optim.minimizer(res)
            log_term2 = log((comp.fluid.d_Hyd/2+comp.membrane.thick)/(comp.fluid.d_Hyd/2))
            comp.J_perm = comp.membrane.D / (comp.fluid.d_Hyd/2 * log_term2) *
                          (comp.membrane.K_S * (cw / comp.fluid.Solubility - p_out^0.5))
            return cw
        else
            function eq_full(cv)
                c_wl, c_ws = cv
                J_mt   = comp.fluid.k_t * (c - c_wl)
                log_t  = log((comp.fluid.d_Hyd/2 + comp.membrane.thick)/(comp.fluid.d_Hyd/2))
                J_d    = comp.membrane.k_d * (c_wl / comp.membrane.K_S) -
                         comp.membrane.k_d * comp.membrane.K_S * c_ws^2
                J_diff = comp.membrane.D / (comp.fluid.d_Hyd/2 * log_t) *
                         (comp.membrane.K_S * c_ws - p_out^0.5)
                return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
            end
            ub = c * (1 + 1e-4)
            res = optimize(eq_full, [0.0, 0.0], [ub, ub], [2c/3, c/3], Fminbox(NelderMead()),
                           Optim.Options(g_tol=1e-8))
            cwl = Optim.minimizer(res)[1]
            comp.J_perm = comp.fluid.k_t * (c - cwl)
            return cwl
        end
    end
    return c_guess
end

# ---------------------------------------------------------------------------
# Numerical efficiency (1-D integration along pipe length)
# ---------------------------------------------------------------------------

"""
    get_efficiency!(comp; plotvar=false, c_guess=nothing, p_out=1e-15)

Numerically integrate the concentration profile along the pipe length
to obtain extraction efficiency.
"""
function get_efficiency!(comp::Component;
                         c_guess::Union{Float64,Nothing}=nothing,
                         p_out::Float64=1e-15)
    if comp.c_in == 0
        comp.c_out = 0.0; comp.eff = 0.0; return
    end
    L_vec = range(0, comp.geometry.L, length=100)
    dl    = step(L_vec)
    c_vec = zeros(Float64, length(L_vec))
    f_H2  = comp.fluid.MS ? 0.5 : 1.0
    cg    = c_guess !== nothing ? c_guess : Float64(comp.c_in)

    for i in eachindex(L_vec)
        if i == 1
            c_vec[1] = Float64(comp.c_in)
            cg = get_flux!(comp, c_vec[1]; c_guess=cg, p_out=p_out)
        else
            c_vec[i] = c_vec[i-1] + f_H2 * comp.J_perm * comp.fluid.d_Hyd * π * dl^2 /
                       comp.fluid.U0 / (π * comp.fluid.d_Hyd^2 / 4 * dl)
            cg = get_flux!(comp, c_vec[i]; c_guess=cg, p_out=p_out)
        end
    end
    comp.eff = (comp.c_in - c_vec[end]) / comp.c_in
end

# ---------------------------------------------------------------------------
# Outlet concentration accounting for recirculation
# ---------------------------------------------------------------------------

"""
    outlet_c_comp!(comp) -> Float64

Compute outlet concentration including recirculation effects.
See Python docstring for physics description of modes.
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
# Global HX coefficient
# ---------------------------------------------------------------------------

"""
    get_global_HX_coeff!(comp; R_conv_sec=0.0)

Compute the overall heat transfer coefficient U for a heat-exchanger component.
Stores result in `comp.U` and `comp.fluid.h_coeff`.
"""
function get_global_HX_coeff!(comp::Component; R_conv_sec::Float64=0.0)
    R_cond = log((comp.fluid.d_Hyd + comp.membrane.thick) / comp.fluid.d_Hyd) /
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
    R_conv_prim = 1.0 / h_prim
    comp.U = 1.0 / (R_conv_prim + R_cond + R_conv_sec)
end

# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------

"""
    get_solid_inventory!(comp; p_out=0.0, flag_an=false) -> Float64

Membrane tritium inventory [mol]. Uses numerical 2-D integration by default,
or the analytical formula when `flag_an=true`.
"""
function get_solid_inventory!(comp::Component; p_out::Float64=0.0, flag_an::Bool=false)::Float64
    flag_an && return analytical_solid_inventory!(comp; p_out=p_out)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)

    r_in  = comp.fluid.d_Hyd / 2
    r_out = r_in + comp.membrane.thick
    L_max = comp.geometry.L

    function integrand(r, L)
        if !comp.fluid.MS
            K_S  = comp.membrane.K_S; d_H = comp.fluid.d_Hyd
            dimless  = 2 * comp.membrane.D * K_S /
                       (comp.fluid.k_t * comp.fluid.Solubility * d_H *
                        log((d_H + 2*comp.membrane.thick) / d_H))
            dimless2 = 2 * comp.membrane.D * K_S /
                       (comp.fluid.Solubility * d_H *
                        log((d_H + 2*comp.membrane.thick) / d_H))
            L_ch  = -dimless / (1 + dimless) * 4 * comp.fluid.k_t / (comp.fluid.U0 * d_H)
            c_ext = p_out^0.5 * K_S
            c_w   = (comp.c_in / comp.fluid.Solubility * K_S) * exp(L_ch * L) /
                    (dimless2 / comp.fluid.k_t + 1) + c_ext
            return (-(c_w - c_ext) * log(r / r_out) / log(r_out / r_in) + c_ext) * 2π * r
        else
            d_H = comp.fluid.d_Hyd
            tau = 4 * comp.fluid.k_t * L / (comp.fluid.U0 * d_H)
            xi  = (1.0 / comp.c_in / comp.fluid.Solubility *
                   (0.5 * comp.membrane.K_S * comp.membrane.D /
                    (comp.fluid.k_t * d_H * log((d_H + 2*comp.membrane.thick)/d_H)))^2)
            beta     = (1/xi + 1)^0.5 + log((1/xi + 1)^0.5 - 1)
            beta_tau = beta - tau - 1.0
            if beta_tau > log(floatmax(Float64))
                w = beta_tau - log(beta_tau)
            else
                w = real(lambertw(complex(exp(beta_tau)), 0))
            end
            alpha  = (1.0 / comp.fluid.Solubility *
                      (0.5 * comp.membrane.D * comp.membrane.K_S /
                       (comp.fluid.k_t * d_H * log((d_H + 2*comp.membrane.thick)/d_H)))^2)
            c_ext  = p_out^0.5 * comp.membrane.K_S
            c_w_l  = alpha * (w^2 + 2w) + alpha * (2 - 2*((w^2 + 2w) + 1)^0.5) + c_ext
            c_w_l  = max(c_w_l, 1e-17)
            return ((-log(r / r_out) / log(r_out / r_in) *
                     ((alpha / comp.fluid.Solubility)^0.5 * w * comp.membrane.K_S - c_ext) +
                     c_ext) * 2π * r)
        end
    end

    result, _ = quadgk(r -> quadgk(L -> integrand(r, L), 0.0, L_max)[1], r_in, r_out)
    comp.membrane.inv = result * comp.geometry.n_pipes
    isnan(comp.membrane.inv) && println("Error: Inventory calculation failed")
    return comp.membrane.inv
end

"""
    analytical_solid_inventory!(comp; p_out=0.0) -> Float64

Analytical membrane inventory formula.
"""
function analytical_solid_inventory!(comp::Component; p_out::Float64=0.0)::Float64
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    circ(r) = π * r^2
    if !comp.fluid.MS
        d_H = comp.fluid.d_Hyd; K_S = comp.membrane.K_S
        dimless  = 2*comp.membrane.D*K_S /
                   (comp.fluid.k_t * comp.fluid.Solubility * d_H *
                    log((d_H + 2*comp.membrane.thick)/d_H))
        dimless2 = 2*comp.membrane.D*K_S /
                   (comp.fluid.Solubility * d_H *
                    log((d_H + 2*comp.membrane.thick)/d_H))
        L_ch = -dimless/(1+dimless) * 4*comp.fluid.k_t/(comp.fluid.U0 * d_H)
        K    = -2π * (comp.c_in / (dimless2/comp.fluid.k_t + 1) / comp.fluid.Solubility * K_S) /
               log((comp.geometry.D/2 + comp.geometry.thick) / (comp.geometry.D/2))
        L_factor = (exp(L_ch * comp.geometry.L) - 1) / L_ch
        K       *= L_factor
        ifun(r)  = 1/4 * r^2 * (2*log(r / (comp.geometry.D/2 + comp.geometry.thick)) - 1)
        c_ext    = p_out^0.5 * K_S
        integral = (K * ifun(comp.geometry.D/2 + comp.geometry.thick) -
                    K * ifun(comp.geometry.D/2)) +
                   comp.geometry.L * c_ext * (circ(comp.geometry.D/2 + comp.geometry.thick) -
                                              circ(comp.geometry.D/2))
        comp.membrane.inv = integral
        return integral
    else
        # MS analytical
        analytical_efficiency!(comp; p_out=p_out)
        d_H = comp.fluid.d_Hyd
        function ms_integral(L_val)
            tau      = 4*comp.fluid.k_t*L_val/(comp.fluid.U0*d_H)
            beta     = (1/comp.xi + 1)^0.5 + log((1/comp.xi + 1)^0.5 - 1)
            beta_tau = beta - tau - 1.0
            if beta_tau > log(floatmax(Float64))
                w = beta_tau - log(beta_tau)
            else
                w = real(lambertw(complex(exp(beta_tau)), 0))
            end
            c_ext = p_out^0.5 * comp.membrane.K_S
            K     = (comp.alpha^0.5 / comp.fluid.Solubility^0.5 *
                     (-beta_tau * (w^2 - w + 1) /
                      (4*comp.fluid.k_t/(comp.fluid.U0*d_H))) * comp.membrane.K_S)
            ifun(r) = 1/4 * r^2 * (2*log(r/(comp.geometry.D/2 + comp.geometry.thick)) - 1)
            return K * ifun(comp.geometry.D/2 + comp.geometry.thick) -
                   K * ifun(comp.geometry.D/2)
        end
        p_out_term = comp.geometry.L * p_out^0.5 * comp.membrane.K_S *
                     (circ(comp.geometry.D/2 + comp.geometry.thick) - circ(comp.geometry.D/2))
        inv = ms_integral(comp.geometry.L) - ms_integral(0.0) + p_out_term
        comp.membrane.inv = inv * comp.geometry.n_pipes
        return inv
    end
end

"""
    get_fluid_inventory!(comp; flag_an=false, p_out=0.0) -> Float64

Fluid-phase tritium inventory [mol].
"""
function get_fluid_inventory!(comp::Component; flag_an::Bool=false, p_out::Float64=0.0)::Float64
    flag_an && return analytical_fluid_inventory!(comp; p_out=p_out)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    r_in = comp.fluid.d_Hyd / 2
    d_H  = comp.fluid.d_Hyd

    function integrand(L)
        if !comp.fluid.MS
            dimless = 2*comp.membrane.D*comp.membrane.K_S /
                      (comp.fluid.k_t * comp.fluid.Solubility * d_H *
                       log((d_H + 2*comp.membrane.thick)/d_H))
            L_ch    = -dimless/(1+dimless) * 4*comp.fluid.k_t/(comp.fluid.U0*d_H)
            c_ext   = p_out^0.5 * comp.fluid.Solubility
            return (comp.c_in - c_ext) * exp(L_ch * L) + c_ext
        else
            tau = 4*comp.fluid.k_t*L/(comp.fluid.U0*d_H)
            xi  = (1.0 / comp.c_in / comp.fluid.Solubility *
                   (0.5 * comp.membrane.K_S * comp.membrane.D /
                    (comp.fluid.k_t * d_H * log((d_H+2*comp.membrane.thick)/d_H)))^2)
            beta     = (1/xi + 1)^0.5 + log((1/xi + 1)^0.5 - 1)
            beta_tau = beta - tau - 1.0
            if beta_tau > log(floatmax(Float64))
                w = beta_tau - log(beta_tau)
            else
                w = real(lambertw(complex(exp(beta_tau)), 0))
            end
            alpha = (1.0 / comp.fluid.Solubility *
                     (0.5 * comp.membrane.D * comp.membrane.K_S /
                      (comp.fluid.k_t * d_H * log((d_H+2*comp.membrane.thick)/d_H)))^2)
            return alpha * (w^2 + 2w)
        end
    end

    result, _ = quadgk(integrand, 0.0, comp.geometry.L)
    comp.fluid.inv = result * π * r_in^2 * comp.geometry.n_pipes
    return comp.fluid.inv
end

"""
    analytical_fluid_inventory!(comp; p_out=0.0) -> Float64

Analytical fluid inventory (LM only; MS falls back to numerical).
"""
function analytical_fluid_inventory!(comp::Component; p_out::Float64=0.0)::Float64
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid;
                                   turbulator=comp.geometry.turbulator)
    if !comp.fluid.MS
        d_H  = comp.fluid.d_Hyd
        dimless  = 2*comp.membrane.D*comp.membrane.K_S /
                   (comp.fluid.k_t*comp.fluid.Solubility*d_H*
                    log((d_H+2*comp.membrane.thick)/d_H))
        L_ch     = -dimless/(1+dimless)*4*comp.fluid.k_t/(comp.fluid.U0*d_H)
        c_ext    = p_out^0.5 * comp.fluid.Solubility
        K        = (comp.c_in - c_ext) * (exp(L_ch*comp.geometry.L) - 1) / L_ch
        circ(r)  = π * r^2
        integral = K * circ(d_H/2) + c_ext * circ(d_H/2) * comp.geometry.L
        comp.fluid.inv = integral * comp.geometry.n_pipes
        return integral
    else
        println("MS fluid integration is done numerically")
        return get_fluid_inventory!(comp; flag_an=false, p_out=p_out)
    end
end

"""
    get_inventory!(comp; flag_an=true, p_out=0.0)

Compute total inventory = fluid + membrane. Stores in `comp.inv`.
"""
function get_inventory!(comp::Component; flag_an::Bool=true, p_out::Float64=0.0)
    get_solid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    get_fluid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    comp.inv = comp.fluid.inv + comp.membrane.inv
end

end # module
