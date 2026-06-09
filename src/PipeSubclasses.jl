"""
    PipeSubclasses

Mutable structs representing the geometry, fluid, membrane, and turbulator
sub-objects used by the Component type.
"""
module PipeSubclasses

using ..TriomaTypes: TriomaClass, update_attribute!
using ..Correlations
using AtomicAndPhysicalConstants: BOLTZMANN_k

export Geometry, Fluid, Membrane, FluidMaterial, SolidMaterial,
       Turbulator, WireCoil, CustomTurbulator,
       get_fluid_volume, get_solid_volume, get_total_volume,
       set_properties_from_fluid_material!, set_properties_from_solid_material!,
       update_T_prop!, get_kt!,
       k_t_correlation, h_t_correlation

# ── Turbulators ────────────────────────────────────────────────────────────────

"""Base turbulator type."""
mutable struct Turbulator <: TriomaClass
    turbulator_type::Union{String,Nothing}
end
Turbulator(; turbulator_type=nothing) = Turbulator(turbulator_type)

"""Wire-coil turbulator."""
mutable struct WireCoil <: TriomaClass
    turbulator_type::String
    pitch::Union{Float64,Nothing}
end
WireCoil(; pitch=nothing) = WireCoil("WireCoil", pitch)

function k_t_correlation(t::WireCoil; Re, Sc, d_hyd, D)
    Sh = Re > 2030 ? 0.132 * Re^0.72 * Sc^0.37 * (t.pitch / d_hyd)^(-0.372) : 3.66
    return Correlations.get_k_from_Sh(Sh, d_hyd, D)
end

function h_t_correlation(t::WireCoil; Re, Pr, d_hyd, k)
    Nu = Re > 2030 ? 0.132 * Re^0.72 * Pr^0.37 * (t.pitch / d_hyd)^(-0.372) : 3.66
    return Correlations.get_h_from_Nu(Nu, k, d_hyd)
end

"""Custom (user-defined power-law) turbulator."""
mutable struct CustomTurbulator <: TriomaClass
    turbulator_type::String
    a::Union{Float64,Nothing}
    b::Union{Float64,Nothing}
    c::Union{Float64,Nothing}
end
CustomTurbulator(; a=nothing, b=nothing, c=nothing) =
    CustomTurbulator("Custom", a, b, c)

function k_t_correlation(t::CustomTurbulator; Re, Sc, d_hyd, D)
    Sh = Re > 2030 ? t.a * Re^t.b * Sc^t.c : 3.66
    return Correlations.get_k_from_Sh(Sh, d_hyd, D)
end

function h_t_correlation(t::CustomTurbulator; Re, Pr, d_hyd, k)
    Nu = Re > 2030 ? t.a * Re^t.b * Pr^t.c : 3.66
    # Note: the Python source references self.pitch here, which appears to be a
    # copy-paste bug (CustomTurbulator has no pitch). Using d_hyd as characteristic
    # length, consistent with the WireCoil implementation.
    return Correlations.get_h_from_Nu(Nu, k, d_hyd)
end

# ── Geometry ───────────────────────────────────────────────────────────────────

"""
    Geometry

Geometric parameters for a pipe/tube component.

Fields
- `L`        : length [m]
- `D`        : inner diameter [m]
- `thick`    : wall thickness [m]
- `n_pipes`  : number of parallel pipes (default 1)
- `turbulator`: optional turbulator object
"""
mutable struct Geometry <: TriomaClass
    L::Union{Float64,Nothing}
    D::Union{Float64,Nothing}
    thick::Union{Float64,Nothing}
    n_pipes::Float64
    turbulator::Union{WireCoil,CustomTurbulator,Turbulator,Nothing}
end
function Geometry(; L=nothing, D=nothing, thick=nothing, n_pipes=1.0, turbulator=nothing)
    Geometry(L, D, thick, n_pipes, turbulator)
end

"""Fluid volume of a single pipe [m³]."""
get_fluid_volume(g::Geometry) = π * (g.D / 2)^2 * g.L

"""Wall (solid) volume of a single pipe [m³]."""
get_solid_volume(g::Geometry) = π * ((g.D / 2)^2 - (g.D / 2 - g.thick)^2) * g.L

"""Total volume (fluid + wall) of a single pipe [m³]."""
get_total_volume(g::Geometry) = get_fluid_volume(g) + get_solid_volume(g)

# ── FluidMaterial / SolidMaterial (simple property containers) ─────────────────

"""Physical properties of a fluid material at a given temperature."""
mutable struct FluidMaterial <: TriomaClass
    T::Union{Float64,Nothing}
    D::Union{Float64,Nothing}          # tritium diffusivity [m²/s]
    Solubility::Union{Float64,Nothing} # Henry or Sievert constant
    MS::Union{Bool,Nothing}            # true = molten salt, false = liquid metal
    mu::Union{Float64,Nothing}         # dynamic viscosity [Pa·s]
    rho::Union{Float64,Nothing}        # density [kg/m³]
    k::Union{Float64,Nothing}          # thermal conductivity [W/m/K]
    cp::Union{Float64,Nothing}         # specific heat [J/kg/K]
end
function FluidMaterial(; T=nothing, D=nothing, Solubility=nothing,
                         MS=nothing, mu=nothing, rho=nothing, k=nothing, cp=nothing)
    FluidMaterial(T, D, Solubility, MS, mu, rho, k, cp)
end

"""Physical properties of a solid (membrane) material at a given temperature."""
mutable struct SolidMaterial <: TriomaClass
    T::Union{Float64,Nothing}
    D::Union{Float64,Nothing}    # tritium diffusivity [m²/s]
    K_S::Union{Float64,Nothing}  # Sievert constant
    k::Union{Float64,Nothing}    # thermal conductivity [W/m/K]
end
SolidMaterial(; T=nothing, D=nothing, K_S=nothing, k=nothing) =
    SolidMaterial(T, D, K_S, k)

# ── Fluid ──────────────────────────────────────────────────────────────────────

const K_B_EV = BOLTZMANN_k   # Boltzmann constant [eV/K]

"""
    Fluid

Tritium-bearing fluid properties for transport analysis.

Temperature-dependent diffusivity and solubility are computed automatically
if Arrhenius parameters (`D_0`/`E_d` or `Solubility_0`/`E_s`) are provided.
"""
mutable struct Fluid <: TriomaClass
    T::Union{Float64,Nothing}
    MS::Bool
    D_0::Union{Float64,Nothing}
    E_d::Union{Float64,Nothing}
    Solubility_0::Union{Float64,Nothing}
    E_s::Union{Float64,Nothing}
    D::Union{Float64,Nothing}
    Solubility::Union{Float64,Nothing}
    k_t::Union{Float64,Nothing}
    d_Hyd::Union{Float64,Nothing}
    mu::Union{Float64,Nothing}
    rho::Union{Float64,Nothing}
    recirculation::Float64
    U0::Union{Float64,Nothing}
    k::Union{Float64,Nothing}
    cp::Union{Float64,Nothing}
    inv::Union{Float64,Nothing}
    V::Union{Float64,Nothing}
    h_coeff::Union{Float64,Nothing}   # heat transfer coefficient [W/m²/K], set by HX calcs
end

function Fluid(;
    T=nothing, D=nothing, D_0=nothing, E_d=nothing,
    Solubility=nothing, Solubility_0=nothing, E_s=nothing,
    MS=true, d_Hyd=nothing, k_t=nothing, mu=nothing, rho=nothing,
    U0=nothing, k=nothing, cp=nothing, inv=nothing, recirculation=0.0, V=nothing
)
    D_val = (D_0 !== nothing && E_d !== nothing && T !== nothing) ?
        D_0 * exp(-E_d / (K_B_EV * T)) : D
    Sol_val = (Solubility_0 !== nothing && E_s !== nothing && T !== nothing) ?
        Solubility_0 * exp(-E_s / (K_B_EV * T)) : Solubility
    U0_eff = (recirculation == 0.0 || U0 === nothing) ? U0 : U0 * (1 + recirculation)
    Fluid(T, MS, D_0, E_d, Solubility_0, E_s, D_val, Sol_val,
          k_t, d_Hyd, mu, rho, recirculation, U0_eff, k, cp, inv, V, nothing)
end

"""Copy fluid thermo-physical properties from a `FluidMaterial` object."""
function set_properties_from_fluid_material!(fluid::Fluid, mat::FluidMaterial)
    fluid.T         = mat.T
    fluid.D         = mat.D
    fluid.Solubility = mat.Solubility
    fluid.mu        = mat.mu
    fluid.rho       = mat.rho
    fluid.cp        = mat.cp
    fluid.k         = mat.k
end

"""Re-evaluate temperature-dependent D and Solubility using Arrhenius expressions."""
function update_T_prop!(fluid::Fluid)
    if fluid.D_0 !== nothing && fluid.E_d !== nothing && fluid.T !== nothing
        fluid.D = fluid.D_0 * exp(-fluid.E_d / (K_B_EV * fluid.T))
    end
    if fluid.Solubility_0 !== nothing && fluid.E_s !== nothing && fluid.T !== nothing
        fluid.Solubility = fluid.Solubility_0 * exp(-fluid.E_s / (K_B_EV * fluid.T))
    end
end

"""
    get_kt!(fluid::Fluid; turbulator=nothing)

Compute the mass-transfer coefficient `k_t` from the fluid's hydraulic diameter,
velocity, and thermophysical properties. Uses the Getthem correlation for turbulent
flow and Sh = 3.66 for laminar flow, or the turbulator's own correlation if one is
provided.
"""
function get_kt!(fluid::Fluid; turbulator=nothing)
    if fluid.d_Hyd === nothing
        @warn "Hydraulic diameter is not defined"
        return
    end
    if fluid.k_t !== nothing
        println("k_t is already defined")
        return
    end
    Re_val = Correlations.Re(fluid.rho, fluid.U0, fluid.d_Hyd, fluid.mu)
    Sc_val = Correlations.Schmidt(fluid.D, fluid.mu, fluid.rho)

    if turbulator === nothing
        if Re_val > 2030
            Sh = 0.0096 * Re_val^0.913 * Sc_val^0.346
        else
            println("Re = $(Re_val) indicates laminar flow")
            Sh = 3.66
        end
        fluid.k_t = Correlations.get_k_from_Sh(Sh, fluid.d_Hyd, fluid.D)
    else
        if turbulator isa WireCoil || turbulator isa CustomTurbulator
            fluid.k_t = k_t_correlation(turbulator; Re=Re_val, Sc=Sc_val,
                                         d_hyd=fluid.d_Hyd, D=fluid.D)
        else
            error("Turbulator type '$(turbulator.turbulator_type)' not implemented")
        end
    end
end

# ── Membrane ───────────────────────────────────────────────────────────────────

"""
    Membrane

Metallic membrane properties for hydrogen-isotope transport.

Temperature-dependent `D` and `K_S` are evaluated from Arrhenius parameters
(`D_0`/`E_d` and `K_S_0`/`E_S`) if provided.
"""
mutable struct Membrane <: TriomaClass
    T::Union{Float64,Nothing}
    D::Union{Float64,Nothing}
    thick::Union{Float64,Nothing}
    K_S::Union{Float64,Nothing}
    k_d::Union{Float64,Nothing}   # dissociation rate constant
    k_r::Union{Float64,Nothing}   # recombination rate constant
    k::Union{Float64,Nothing}     # thermal conductivity [W/m/K]
    D_0::Union{Float64,Nothing}
    E_d::Union{Float64,Nothing}
    K_S_0::Union{Float64,Nothing}
    E_S::Union{Float64,Nothing}
    inv::Union{Float64,Nothing}
    V::Union{Float64,Nothing}
end

function Membrane(;
    T=nothing, D=nothing, thick=nothing, K_S=nothing,
    k_d=nothing, k_r=nothing, k=nothing,
    D_0=nothing, E_d=nothing, K_S_0=nothing, E_S=nothing,
    inv=nothing, V=nothing
)
    D_val  = (D_0 !== nothing && E_d !== nothing && T !== nothing) ?
        D_0 * exp(-E_d / (K_B_EV * T)) : D
    KS_val = (K_S_0 !== nothing && E_S !== nothing && T !== nothing) ?
        K_S_0 * exp(-E_S / (K_B_EV * T)) : K_S
    Membrane(T, D_val, thick, KS_val, k_d, k_r, k, D_0, E_d, K_S_0, E_S, inv, V)
end

"""Copy membrane properties from a `SolidMaterial` object."""
function set_properties_from_solid_material!(mem::Membrane, mat::SolidMaterial)
    mem.T   = mat.T
    mem.D   = mat.D
    mem.K_S = mat.K_S
end

"""Re-evaluate temperature-dependent D and K_S using Arrhenius expressions."""
function update_T_prop!(mem::Membrane)
    if mem.D_0 !== nothing && mem.E_d !== nothing && mem.T !== nothing
        mem.D = mem.D_0 * exp(-mem.E_d / (K_B_EV * mem.T))
    end
    if mem.K_S_0 !== nothing && mem.E_S !== nothing && mem.T !== nothing
        mem.K_S = mem.K_S_0 * exp(-mem.E_S / (K_B_EV * mem.T))
    end
end

end # module PipeSubclasses
