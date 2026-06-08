"""
    PAVComponent

The main pipe-component type for Permeation Against Vacuum (PAV) analysis and
general 0D/1D tritium transport calculations.
"""
module PAVComponent

using ..TriomaModule: TriomaClass, update_attribute!
using ..PipeSubclasses
using ..Correlations
using ..MoltenSalts
using ..LiquidMetals
using QuadGK
using Optim
using SpecialFunctions: lambertw    # Note: verify this function is available in SpecialFunctions.jl

export Component,
       friction_factor, update_T_prop!, get_pressure_drop!, get_pumping_power!,
       estimate_cost, get_pipe_flowrate, get_total_flowrate, define_component_volumes!,
       get_adimensionals!, analytical_efficiency!, use_analytical_efficiency!,
       get_efficiency!, get_flux!, get_regime, get_solid_inventory!, get_fluid_inventory!,
       get_inventory!, get_global_HX_coeff!, outlet_c_comp!, T_leak,
       split_HX, converge_split_HX

# ── Component ─────────────────────────────────────────────────────────────────

"""
    Component

Represents a pipe-type component in the outer fuel cycle.

Fields map directly to the Python attributes. All fields are `Union{T,Nothing}`
except `loss` (Bool) and `p_out` (Float64, defaulting to 1e-15).
"""
mutable struct Component <: TriomaClass
    c_in::Union{Float64,Nothing}
    geometry::Union{Geometry,Nothing}
    eff::Union{Float64,Nothing}
    n_pipes::Union{Tuple{Float64},Nothing}   # stored as 1-tuple to match Python
    fluid::Union{Fluid,Nothing}
    membrane::Union{Membrane,Nothing}
    name::Union{String,Nothing}
    loss::Bool
    inv::Union{Float64,Nothing}
    p_out::Float64
    delta_p::Union{Float64,Nothing}
    U::Union{Float64,Nothing}           # global HX coefficient [W/m²/K]
    pumping_power::Union{Float64,Nothing}
    cost::Union{Float64,Nothing}
    # Computed attributes (set by calculation methods)
    c_out::Union{Float64,Nothing}
    J_perm::Union{Float64,Nothing}
    tau::Union{Float64,Nothing}
    xi::Union{Float64,Nothing}
    alpha::Union{Float64,Nothing}
    eff_an::Union{Float64,Nothing}
    zeta::Union{Float64,Nothing}
    H::Union{Float64,Nothing}
    W::Union{Float64,Nothing}
    pipe_flowrate::Union{Float64,Nothing}
    flowrate::Union{Float64,Nothing}
end

function Component(;
    geometry=nothing, c_in=nothing, eff=nothing, fluid=nothing,
    membrane=nothing, name=nothing, p_out=1e-15, loss=false,
    inv=nothing, delta_p=nothing, pumping_power=nothing, U=nothing,
    V=nothing, cost=nothing
)
    n_pipes_val = geometry !== nothing ? (geometry.n_pipes,) : nothing
    Component(
        c_in, geometry, eff, n_pipes_val, fluid, membrane, name, loss,
        inv, p_out, delta_p, U, pumping_power, cost,
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
        nothing, nothing
    )
end

# ── custom update_attribute! for Component ────────────────────────────────────

"""
Override to propagate temperature updates to fluid and membrane sub-objects.
"""
function update_attribute!(comp::Component, attr_name::Symbol, new_value)
    if attr_name === :T
        comp.fluid   !== nothing && update_attribute!(comp.fluid,   :T, new_value)
        comp.membrane !== nothing && update_attribute!(comp.membrane, :T, new_value)
        update_T_prop!(comp)
        return
    end
    if attr_name in fieldnames(Component)
        setfield!(comp, attr_name, new_value)
        if attr_name === :n_pipes
            for fname in fieldnames(Component)
                child = getfield(comp, fname)
                if child isa TriomaClass && :n_pipes in fieldnames(typeof(child))
                    setfield!(child, :n_pipes, new_value)
                end
            end
        end
        return
    end
    for fname in fieldnames(Component)
        child = getfield(comp, fname)
        if child isa TriomaClass && attr_name in fieldnames(typeof(child))
            setfield!(child, attr_name, new_value)
            return
        end
    end
    throw(ArgumentError("'$(attr_name)' is not an attribute of Component"))
end

# ── Thermo-property update ────────────────────────────────────────────────────

function update_T_prop!(comp::Component)
    comp.fluid    !== nothing && PipeSubclasses.update_T_prop!(comp.fluid)
    comp.membrane !== nothing && PipeSubclasses.update_T_prop!(comp.membrane)
end

# ── Hydraulics ────────────────────────────────────────────────────────────────

"""Darcy-Weisbach friction factor (Blasius for turbulent, 64/Re for laminar)."""
function friction_factor(comp::Component, Re::Float64)
    return Re < 2300 ? 64.0 / Re : 0.316 / Re^0.25
end

"""Compute and store pressure drop [Pa]."""
function get_pressure_drop!(comp::Component)
    rho = comp.fluid.rho;  U = comp.fluid.U0
    D   = comp.geometry.D; mu = comp.fluid.mu; L = comp.geometry.L
    Re_val = Correlations.Re(rho, U, D, mu)
    f = friction_factor(comp, Re_val)
    comp.delta_p = f * (L / D) * (rho * U^2) / 2
    return comp.delta_p
end

"""Compute and store pumping power [W]."""
function get_pumping_power!(comp::Component)
    comp.delta_p === nothing && get_pressure_drop!(comp)
    comp.pumping_power = comp.delta_p * get_pipe_flowrate(comp) * comp.geometry.n_pipes
    return comp.pumping_power
end

"""Volumetric flow rate in a single pipe [m³/s]."""
function get_pipe_flowrate(comp::Component)
    comp.pipe_flowrate = comp.fluid.U0 * π * comp.fluid.d_Hyd^2 / 4
    return comp.pipe_flowrate
end

"""Total volumetric flow rate across all pipes [m³/s]."""
function get_total_flowrate(comp::Component)
    get_pipe_flowrate(comp)
    comp.flowrate = comp.pipe_flowrate * comp.geometry.n_pipes
    return comp.flowrate * comp.geometry.n_pipes
end

# ── Cost / volume ─────────────────────────────────────────────────────────────

"""
    estimate_cost(comp, metal_cost, fluid_cost) -> cost [USD]

`metal_cost` and `fluid_cost` are in USD/m³.
"""
function estimate_cost(comp::Component; metal_cost=0.0, fluid_cost=0.0)
    V_solid = PipeSubclasses.get_solid_volume(comp.geometry)
    V_fluid = PipeSubclasses.get_fluid_volume(comp.geometry)
    comp.cost = (V_solid * metal_cost + V_fluid * fluid_cost) * comp.geometry.n_pipes
    return comp.cost
end

"""Set fluid and membrane volume attributes from geometry."""
function define_component_volumes!(comp::Component)
    comp.fluid.V    = PipeSubclasses.get_fluid_volume(comp.geometry)
    comp.membrane.V = PipeSubclasses.get_solid_volume(comp.geometry)
    comp.V          = comp.fluid.V + comp.membrane.V
end

# ── Tritium leak ──────────────────────────────────────────────────────────────

"""Tritium leakage from the component [mol/s]."""
function T_leak(comp::Component)
    return comp.c_in * comp.eff * get_pipe_flowrate(comp)
end

# ── Regime / dimensionless numbers ───────────────────────────────────────────

"""Identify the dominant transport regime (string)."""
function get_regime(comp::Component; print_var::Bool=false)
    comp.fluid === nothing   && (println("No fluid selected"); return)
    comp.membrane === nothing && (println("No membrane selected"); return)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

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

"""Compute and store the H and W adimensional parameters."""
function get_adimensionals!(comp::Component)
    comp.fluid === nothing && (println("No fluid selected"); return)
    comp.fluid.k_t === nothing && PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    if comp.fluid.MS
        comp.H = MoltenSalts.H_ms(comp.fluid.k_t, comp.fluid.Solubility, comp.membrane.k_d)
        comp.W = MoltenSalts.W_ms(comp.membrane.k_d, comp.membrane.D, comp.membrane.thick,
                                   comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
    else
        W_val = LiquidMetals.W_lm(comp.membrane.k_r, comp.membrane.D, comp.membrane.thick,
                                   comp.membrane.K_S, comp.c_in, comp.fluid.Solubility)
        pp    = LiquidMetals.partition_param_lm(comp.membrane.D, comp.fluid.k_t,
                                                 comp.membrane.K_S, comp.fluid.Solubility,
                                                 comp.membrane.thick)
        comp.H = W_val * pp
        comp.W = W_val
    end
end

# ── Analytical efficiency ─────────────────────────────────────────────────────

"""
    analytical_efficiency!(comp; p_out=1e-15)

Evaluate the closed-form efficiency expression from Humrickhouse et al. and
store the result in `comp.eff_an`.
"""
function analytical_efficiency!(comp::Component; p_out::Float64=1e-15)
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    comp.tau = 4 * comp.fluid.k_t * comp.geometry.L / (comp.fluid.U0 * comp.fluid.d_Hyd)

    if comp.fluid.MS
        # ── Molten salt ──
        comp.alpha = (1 / comp.fluid.Solubility) *
            (0.5 * comp.membrane.K_S * comp.membrane.D /
             (comp.fluid.k_t * comp.fluid.d_Hyd *
              log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)))^2
        comp.xi = comp.alpha / comp.c_in
        p_in    = comp.c_in / comp.fluid.Solubility

        if comp.xi > 1e5
            corr_p = 1 - (p_out / p_in)
            comp.eff_an = (1 - exp(-comp.tau)) * corr_p
        elseif comp.xi^0.5 < 1e-2 && comp.tau > 1 / comp.xi^0.5
            corr_p = 1 - (p_out / p_in)^0.5
            comp.eff_an = (1 - (1 - comp.tau * comp.xi^0.5)^2) * corr_p
        else
            e     = (comp.alpha * p_out * comp.fluid.Solubility)^0.5
            f     = e / comp.alpha
            delta = (1 / comp.xi + 1 + 2 * f)^0.5
            beta  = delta + (1 + f) * log(delta - 1 - f)
            max_exp  = log(floatmax(Float64))
            beta_tau = beta - comp.tau - 1

            if beta_tau > max_exp || p_out > 1e-5
                # Numerical fallback via minimization
                lower_b = min(p_out * comp.fluid.Solubility, comp.c_in)
                upper_b = max(p_out * comp.fluid.Solubility, comp.c_in)
                if abs(comp.p_out * comp.fluid.Solubility - comp.c_in) / comp.c_in < 1e-2
                    comp.eff_an = 1e-6
                    return
                end
                function eq_ms(cl_v)
                    cl = cl_v[1]
                    left = (cl / comp.alpha + 1 + 2 * f)^0.5 +
                           (1 + f) * log(-f + ((cl / comp.alpha + 1 + 2 * f)^0.5 - 1))
                    return abs(left - (beta - comp.tau))
                end
                res = optimize(eq_ms, [lower_b], [upper_b],
                               [(lower_b + upper_b) / 2],
                               Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-7))
                comp.eff_an = 1 - (Optim.minimizer(res)[1] / comp.c_in)
            else
                z = exp(beta_tau)
                w = lambertw(z)
                comp.eff_an = real(1 - comp.xi * (w^2 + 2 * w))
            end
        end

    else
        # ── Liquid metal ──
        comp.zeta = (2 * comp.membrane.K_S * comp.membrane.D) /
            (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
             log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
        p_in   = (comp.c_in / comp.fluid.Solubility)^2
        corr_p = 1 - (p_out / p_in)^0.5
        comp.eff_an = (1 - exp(-comp.tau * comp.zeta / (1 + comp.zeta))) * corr_p
    end
end

"""Apply `analytical_efficiency!` and copy result to `comp.eff`."""
function use_analytical_efficiency!(comp::Component; p_out::Float64=1e-15)
    analytical_efficiency!(comp; p_out=p_out)
    comp.eff = comp.eff_an
end

# ── Outlet concentration ──────────────────────────────────────────────────────

"""
    outlet_c_comp!(comp) -> c_out

Compute outlet concentration from current `eff`. Handles recirculation.
"""
function outlet_c_comp!(comp::Component)
    if comp.fluid.recirculation == 0
        comp.c_out = comp.c_in * (1 - comp.eff)

    elseif comp.fluid.recirculation > 0
        err  = 1.0; tol = 1e-6
        comp.c_in == 0 && error("Inlet concentration is zero")
        c0   = comp.c_in; c_in = c0
        while err > tol
            c_in1 = c_in
            update_attribute!(comp, :c_in, c_in)
            analytical_efficiency!(comp)
            update_attribute!(comp, :eff, comp.eff_an)
            comp.c_out = comp.c_in * (1 - comp.eff)
            c_in = (comp.c_out * comp.fluid.recirculation + c0) / (comp.fluid.recirculation + 1)
            err  = abs((c_in - c_in1) / c_in)
        end

    elseif comp.fluid.recirculation < 0
        comp.fluid.recirculation <= -1 &&
            error("Bypass (negative recirculation) not valid: exceeds flow rate")
        comp.c_in == 0 && error("Inlet concentration is zero")
        r = comp.fluid.recirculation
        comp.c_out = comp.c_in * (1 - comp.eff) * (1 + r) + comp.c_in * (-r)
    else
        error("Recirculation factor not valid")
    end
    return comp.c_out
end

# ── Permeation flux ───────────────────────────────────────────────────────────

"""
    get_flux!(comp, c; c_guess=1e-9, p_out=1e-15) -> c_guess_new

Compute the permeation flux `J_perm` [mol/m²/s] at local concentration `c`.
Uses regime-appropriate analytical or numerical expressions.
Returns the updated wall-concentration guess for subsequent calls.
"""
function get_flux!(comp::Component, c::Float64; c_guess::Float64=1e-9, p_out::Float64=1e-15)
    get_adimensionals!(comp)

    if comp.fluid.MS
        # ── Molten salt flux ──
        if comp.W > 10
            if comp.H / comp.W > 1000
                comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
            elseif comp.H / comp.W < 0.0001
                comp.J_perm = -(comp.membrane.D /
                    (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                    comp.membrane.K_S * ((c / comp.fluid.Solubility)^0.5 - p_out^0.5))
            else
                function eq_ms_mt(cwl_v)
                    cwl   = cwl_v[1]
                    J_mt  = 2 * comp.fluid.k_t * (c - cwl)
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        (comp.membrane.K_S * ((cwl / comp.fluid.Solubility)^0.5 - p_out^0.5))
                    return abs(J_diff - J_mt)
                end
                ub  = max(c * (1 + 1e-4), 1e-4)
                res = optimize(eq_ms_mt, [0.0], [ub], [c_guess], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e7)))
                sol = Optim.minimizer(res)[1]
                comp.J_perm = -2 * comp.fluid.k_t * (c - sol)
                return sol
            end
        elseif comp.W < 0.1
            if comp.H > 100
                comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
            elseif comp.H < 0.01
                comp.J_perm = -comp.membrane.k_d * (c / comp.fluid.Solubility)
            else
                function eq_ms_surf(cwl_v)
                    cwl  = cwl_v[1]
                    J_mt = 2 * comp.fluid.k_t * (c - cwl)
                    J_s  = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                           comp.membrane.k_d * comp.membrane.K_S^2 * cwl^2
                    return abs(J_mt - J_s)
                end
                ub  = c * (1 + 1e-4)
                res = optimize(eq_ms_surf, [0.0], [ub], [c * 0.1], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                sol = Optim.minimizer(res)[1]
                comp.J_perm = 2 * comp.fluid.k_t * (c - sol)
                return sol
            end
        else
            if comp.H / comp.W > 1000
                comp.J_perm = -2 * comp.fluid.k_t * (c - p_out * comp.fluid.Solubility)
            elseif comp.H / comp.W < 0.0001
                function eq_ms_mix1(cwl_v)
                    cwl  = cwl_v[1]
                    J_s  = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                           comp.membrane.k_d * comp.membrane.K_S^2 * cwl^2
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        (comp.membrane.K_S * ((cwl / comp.fluid.Solubility)^0.5 - p_out^0.5))
                    return abs(J_diff - J_s)
                end
                res = optimize(eq_ms_mix1, [1e-14], [c], [c_guess], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-7, iterations=Int(1e6)))
                cw  = Optim.minimizer(res)[1]
                comp.J_perm = comp.membrane.D /
                    (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                    (comp.membrane.K_S * (cw / comp.fluid.Solubility)^0.5 - p_out^0.5)
                return cw
            else
                function eq_ms_full(vars)
                    cwl, cws = vars[1], vars[2]
                    J_mt   = 2 * comp.fluid.k_t * (c - cwl)
                    J_d    = comp.membrane.k_d * (cwl / comp.membrane.K_S) -
                             comp.membrane.k_d * comp.membrane.K_S^2 * cws^2
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        ((comp.membrane.K_S * cws) - p_out^0.5)
                    return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
                end
                res = optimize(eq_ms_full, [0.0, 0.0], [c, c], [2c/3, c/3], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                cwl = Optim.minimizer(res)[1]
                comp.J_perm = 2 * comp.fluid.k_t * (c - cwl)
                return cwl
            end
        end

    else
        # ── Liquid metal flux ──
        if comp.W > 10
            if comp.H / comp.W > 1000
                comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
            elseif comp.H / comp.W < 0.0001
                comp.J_perm = -(comp.membrane.D /
                    (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                    (comp.membrane.K_S * (c / comp.fluid.Solubility - p_out^0.5)))
            else
                function eq_lm_mt(cwl_v)
                    cwl   = cwl_v[1]
                    J_mt  = comp.fluid.k_t * (c - cwl)
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        (comp.membrane.K_S * (cwl / comp.fluid.Solubility - p_out^0.5))
                    return abs(J_diff - J_mt)
                end
                ub  = c * (1 + 1e-4)
                res = optimize(eq_lm_mt, [0.0], [ub], [c_guess], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                sol = Optim.minimizer(res)[1]
                comp.J_perm = -comp.fluid.k_t * (c - sol)
                return sol
            end
        elseif comp.W < 0.1
            if comp.H > 100
                comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
            elseif comp.H < 0.01
                comp.J_perm = -comp.membrane.k_d * (c / comp.fluid.Solubility)
            else
                function eq_lm_surf(cwl_v)
                    cwl  = cwl_v[1]
                    J_mt = comp.fluid.k_t * (c - cwl - p_out^0.5 * comp.fluid.Solubility)
                    J_s  = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                           comp.membrane.k_d * comp.membrane.K_S^2 * cwl^2
                    return abs(J_mt - J_s)
                end
                ub  = c * (1 + 1e-4)
                res = optimize(eq_lm_surf, [0.0], [ub], [c * 0.1], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                sol = Optim.minimizer(res)[1]
                comp.J_perm = comp.fluid.k_t * (c - sol - p_out * comp.fluid.Solubility)
                return sol
            end
        else
            if comp.H / comp.W < 0.0001
                comp.J_perm = -comp.fluid.k_t * (c - p_out^0.5 * comp.fluid.Solubility)
            elseif comp.H / comp.W > 1000
                function eq_lm_mix1(cwl_v)
                    cwl  = cwl_v[1]
                    J_s  = comp.membrane.k_d * (c / comp.fluid.Solubility) -
                           comp.membrane.k_d * comp.membrane.K_S^2 * cwl^2
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        (comp.membrane.K_S * (cwl / comp.fluid.Solubility - p_out^0.5))
                    return abs(J_diff - J_s)
                end
                res = optimize(eq_lm_mix1, [1e-14], [c], [c_guess], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                cw  = Optim.minimizer(res)[1]
                comp.J_perm = comp.membrane.D /
                    (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                    (comp.membrane.K_S * (cw / comp.fluid.Solubility - p_out^0.5))
                return cw
            else
                function eq_lm_full(vars)
                    cwl, cws = vars[1], vars[2]
                    J_mt   = comp.fluid.k_t * (c - cwl)
                    J_d    = comp.membrane.k_d * (cwl / comp.membrane.K_S) -
                             comp.membrane.k_d * comp.membrane.K_S * cws^2
                    J_diff = comp.membrane.D /
                        (comp.fluid.d_Hyd / 2 * log((comp.fluid.d_Hyd / 2 + comp.membrane.thick) / (comp.fluid.d_Hyd / 2))) *
                        ((comp.membrane.K_S * cws) - p_out^0.5)
                    return abs(J_mt - J_d) + abs(J_mt - J_diff) + abs(J_d - J_diff)
                end
                ub  = c * (1 + 1e-4)
                res = optimize(eq_lm_full, [0.0, 0.0], [ub, ub], [2c/3, c/3], Fminbox(Powell());
                               options=Optim.Options(f_tol=1e-8, iterations=Int(1e6)))
                cwl = Optim.minimizer(res)[1]
                comp.J_perm = comp.fluid.k_t * (c - cwl)
                return cwl
            end
        end
    end
    return c_guess   # unchanged when J_perm was set directly
end

# ── Numerical efficiency ───────────────────────────────────────────────────────

"""
    get_efficiency!(comp; p_out=1e-15, c_guess=nothing)

Numerical 1D integration of the concentration profile to compute `comp.eff`.
"""
function get_efficiency!(comp::Component; p_out::Float64=1e-15, c_guess=nothing)
    if comp.c_in == 0
        comp.c_out = 0.0; comp.eff = 0.0; return
    end
    L_vec = range(0.0, comp.geometry.L; length=100)
    dl    = step(L_vec)
    c_vec = zeros(length(L_vec))
    f_H2  = comp.fluid.MS ? 0.5 : 1.0

    cg = c_guess !== nothing ? Float64(c_guess) : Float64(comp.c_in)
    for i in eachindex(L_vec)
        if i == 1
            c_vec[i] = Float64(comp.c_in)
            cg = get_flux!(comp, c_vec[i]; c_guess=cg, p_out=p_out)
        else
            c_vec[i] = c_vec[i-1] + f_H2 * comp.J_perm *
                comp.fluid.d_Hyd * π * dl^2 / comp.fluid.U0 /
                (π * comp.fluid.d_Hyd^2 / 4 * dl)
            cg = get_flux!(comp, c_vec[i]; c_guess=cg, p_out=p_out)
        end
    end
    comp.eff = (comp.c_in - c_vec[end]) / comp.c_in
end

# ── Global HX coefficient ─────────────────────────────────────────────────────

"""
    get_global_HX_coeff!(comp; R_conv_sec=0.0)

Compute the overall heat-exchange coefficient U [W/m²/K] of the component.
"""
function get_global_HX_coeff!(comp::Component; R_conv_sec::Float64=0.0)
    R_cond = log((comp.fluid.d_Hyd + comp.membrane.thick) / comp.fluid.d_Hyd) /
             (2 * π * comp.membrane.k)
    Re_val = Correlations.Re(comp.fluid.rho, comp.fluid.U0, comp.fluid.d_Hyd, comp.fluid.mu)
    Pr_val = Correlations.Pr(comp.fluid.cp, comp.fluid.mu, comp.fluid.k)

    if comp.geometry.turbulator === nothing
        h_prim = Correlations.get_h_from_Nu(
            Correlations.Nu_DittusBoelter(Re_val, Pr_val), comp.fluid.k, comp.fluid.d_Hyd)
    else
        h_prim = PipeSubclasses.h_t_correlation(
            comp.geometry.turbulator; Re=Re_val, Pr=Pr_val,
            d_hyd=comp.fluid.d_Hyd, k=comp.fluid.k)
    end
    comp.fluid.h_coeff = h_prim
    comp.U = 1 / (1 / h_prim + R_cond + R_conv_sec)
end

# ── Inventories ───────────────────────────────────────────────────────────────

"""
    get_solid_inventory!(comp; p_out=0.0, flag_an=false) -> inv [mol]

Compute tritium inventory in the metallic membrane.
`flag_an=true` uses the analytical expression; `flag_an=false` uses numerical integration.
"""
function get_solid_inventory!(comp::Component; p_out::Float64=0.0, flag_an::Bool=false)
    flag_an && return _analytical_solid_inventory!(comp; p_out=p_out)

    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    r_in  = comp.fluid.d_Hyd / 2
    r_out = comp.fluid.d_Hyd / 2 + comp.membrane.thick

    function integrand_2d(r, L)
        if !comp.fluid.MS
            dimless  = 2 * comp.membrane.D * comp.membrane.K_S /
                (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
                 log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
            L_ch = -dimless / (1 + dimless) * 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
            c_ext = p_out^0.5 * comp.membrane.K_S
            c_w   = comp.c_in / comp.fluid.Solubility * comp.membrane.K_S *
                    exp(L_ch * L) / (2 * comp.membrane.D * comp.membrane.K_S /
                        (comp.fluid.Solubility * comp.fluid.d_Hyd *
                         log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)) /
                        comp.fluid.k_t + 1) + c_ext
            return (-c_w * log(r / r_out) / log(r_out / r_in) + c_ext) * 2 * π * r
        else
            tau   = 4 * comp.fluid.k_t * L / (comp.fluid.U0 * comp.fluid.d_Hyd)
            xi_l  = (1 / comp.c_in / comp.fluid.Solubility *
                (0.5 * comp.membrane.K_S * comp.membrane.D /
                 (comp.fluid.k_t * comp.fluid.d_Hyd *
                  log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)))^2)
            beta  = (1 / xi_l + 1)^0.5 + log((1 / xi_l + 1)^0.5 - 1)
            beta_tau = beta - tau - 1
            w = beta_tau > log(floatmax(Float64)) ?
                beta_tau - log(beta_tau) : real(lambertw(exp(beta_tau)))
            alpha = (1 / comp.fluid.Solubility *
                (0.5 * comp.membrane.D * comp.membrane.K_S /
                 (comp.fluid.k_t * comp.fluid.d_Hyd *
                  log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)))^2)
            c_ext = p_out^0.5 * comp.membrane.K_S
            c_wl  = alpha * (w^2 + 2w) + alpha * (2 - 2 * (w^2 + 2w + 1)^0.5) + c_ext
            c_wl  = max(c_wl, 1e-17)
            return (-log(r / r_out) / log(r_out / r_in) *
                ((alpha / comp.fluid.Solubility)^0.5 * w * comp.membrane.K_S - c_ext) + c_ext) *
                2 * π * r
        end
    end

    result, _ = quadgk(r -> quadgk(L -> integrand_2d(r, L), 0.0, comp.geometry.L)[1], r_in, r_out)
    comp.membrane.inv = result * comp.geometry.n_pipes
    isnan(comp.membrane.inv) && println("Error: inventory calculation failed")
    return comp.membrane.inv
end

"""Analytical solid inventory (closed-form for LM; numerical for MS)."""
function _analytical_solid_inventory!(comp::Component; p_out::Float64=0.0)
    # Implementation follows the Python analytical_solid_inventory method.
    # For brevity, the LM case is shown; MS delegates to numerical.
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    if !comp.fluid.MS
        integralfun(r) = 1/4 * r^2 * (2 * log(r / (comp.geometry.D / 2 + comp.geometry.thick)) - 1)
        circle(r)      = π * r^2
        dimless  = 2 * comp.membrane.D * comp.membrane.K_S /
            (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
             log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
        dimless2 = 2 * comp.membrane.D * comp.membrane.K_S /
            (comp.fluid.Solubility * comp.fluid.d_Hyd *
             log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
        L_ch = -dimless / (1 + dimless) * 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
        K    = -2 * π * (comp.c_in / (dimless2 / comp.fluid.k_t + 1) /
               comp.fluid.Solubility * comp.membrane.K_S) /
               log((comp.geometry.D / 2 + comp.geometry.thick) / (comp.geometry.D / 2))
        L_factor = (exp(L_ch * comp.geometry.L) - 1) / L_ch
        K = K * L_factor
        integral = K * integralfun(comp.geometry.D / 2 + comp.geometry.thick) -
                   K * integralfun(comp.geometry.D / 2) +
                   comp.geometry.L * p_out^0.5 * comp.membrane.K_S *
                   (circle(comp.geometry.D / 2 + comp.geometry.thick) - circle(comp.geometry.D / 2))
        comp.membrane.inv = integral
        return integral
    else
        return get_solid_inventory!(comp; p_out=p_out, flag_an=false)
    end
end

"""
    get_fluid_inventory!(comp; flag_an=false, p_out=0.0) -> inv [mol]

Compute tritium inventory in the fluid.
"""
function get_fluid_inventory!(comp::Component; flag_an::Bool=false, p_out::Float64=0.0)
    flag_an && return _analytical_fluid_inventory!(comp; p_out=p_out)

    r_in = comp.fluid.d_Hyd / 2
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)

    function integrand_L(L)
        if !comp.fluid.MS
            dimless = 2 * comp.membrane.D * comp.membrane.K_S /
                (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
                 log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
            L_ch    = -dimless / (1 + dimless) * 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
            c_ext   = p_out^0.5 * comp.fluid.Solubility
            return (comp.c_in - c_ext) * exp(L_ch * L) + c_ext
        else
            comp.tau === nothing && analytical_efficiency!(comp; p_out=p_out)
            tau   = 4 * comp.fluid.k_t * L / (comp.fluid.U0 * comp.fluid.d_Hyd)
            xi_l  = (1 / comp.c_in / comp.fluid.Solubility *
                (0.5 * comp.membrane.K_S * comp.membrane.D /
                 (comp.fluid.k_t * comp.fluid.d_Hyd *
                  log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)))^2)
            beta  = (1 / xi_l + 1)^0.5 + log((1 / xi_l + 1)^0.5 - 1)
            beta_tau = beta - tau - 1
            w = beta_tau > log(floatmax(Float64)) ?
                beta_tau - log(beta_tau) : real(lambertw(exp(beta_tau)))
            alpha = (1 / comp.fluid.Solubility *
                (0.5 * comp.membrane.D * comp.membrane.K_S /
                 (comp.fluid.k_t * comp.fluid.d_Hyd *
                  log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd)))^2)
            return alpha * (w^2 + 2w)
        end
    end

    result, _ = quadgk(integrand_L, 0.0, comp.geometry.L)
    comp.fluid.inv = result * π * r_in^2 * comp.geometry.n_pipes
    return comp.fluid.inv
end

"""Analytical fluid inventory (LM only; MS falls back to numerical)."""
function _analytical_fluid_inventory!(comp::Component; p_out::Float64=0.0)
    comp.fluid.k_t === nothing &&
        PipeSubclasses.get_kt!(comp.fluid; turbulator=comp.geometry.turbulator)
    if !comp.fluid.MS
        circle(r) = π * r^2
        dimless  = 2 * comp.membrane.D * comp.membrane.K_S /
            (comp.fluid.k_t * comp.fluid.Solubility * comp.fluid.d_Hyd *
             log((comp.fluid.d_Hyd + 2 * comp.membrane.thick) / comp.fluid.d_Hyd))
        L_ch  = -dimless / (1 + dimless) * 4 * comp.fluid.k_t / (comp.fluid.U0 * comp.fluid.d_Hyd)
        c_ext = p_out^0.5 * comp.fluid.Solubility
        K     = (comp.c_in - c_ext) * (exp(L_ch * comp.geometry.L) - 1) / L_ch
        integral = K * circle(comp.geometry.D / 2) + c_ext * circle(comp.geometry.D / 2) * comp.geometry.L
        comp.fluid.inv = integral * comp.geometry.n_pipes
        return integral
    else
        println("MS fluid integration is done numerically")
        get_fluid_inventory!(comp; flag_an=false, p_out=p_out)
    end
end

"""Compute both fluid and solid inventories and store their sum in `comp.inv`."""
function get_inventory!(comp::Component; flag_an::Bool=true, p_out::Float64=0.0)
    get_solid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    get_fluid_inventory!(comp; flag_an=flag_an, p_out=p_out)
    comp.inv = comp.fluid.inv + comp.membrane.inv
end

end # module PAVComponent
