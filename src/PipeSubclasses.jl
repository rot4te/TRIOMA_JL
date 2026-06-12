"""
    PipeSubclasses

Mutable structs representing the geometry, fluid, membrane, and turbulator
sub-objects used by the Component type.
"""
module PipeSubclasses

using ..TriomaCore: update_attribute!
using ..Correlations
import ..Correlations: hydraulic_diameter
using AtomicAndPhysicalConstants: BOLTZMANN_k

export Geometry, Fluid, Membrane, FluidMaterial, SolidMaterial,
       Turbulator, WireCoil, CustomTurbulator,
       CrossSection, Circular, Rectangular, Annulus, TwistedElliptical,
       flow_area, wetted_perimeter, hydraulic_diameter,
       get_fluid_volume, get_solid_volume, get_total_volume,
       set_properties_from_fluid_material!, set_properties_from_solid_material!,
       update_T_prop!, get_kt!,
       k_t_correlation, h_t_correlation

# ── Turbulators ────────────────────────────────────────────────────────────────

"""Base turbulator type."""
mutable struct Turbulator
    turbulator_type::Union{String,Nothing}
end
Turbulator(; turbulator_type=nothing) = Turbulator(turbulator_type)

"""Wire-coil turbulator."""
mutable struct WireCoil
    turbulator_type::String
    pitch::Union{Float64,Nothing}
end
WireCoil(; pitch=nothing) = WireCoil("WireCoil", pitch)

function k_t_correlation(t::WireCoil;
                         Re::Float64, Sc::Float64,
                         d_hyd::Float64, D::Float64)::Float64
    Sh = Re > 2030 ? 0.132 * Re^0.72 * Sc^0.37 * (t.pitch / d_hyd)^(-0.372) : 3.66
    return Correlations.get_k_from_Sh(Sh, d_hyd, D)
end

function h_t_correlation(t::WireCoil;
                         Re::Float64, Pr::Float64,
                         d_hyd::Float64, k::Float64)::Float64
    Nu = Re > 2030 ? 0.132 * Re^0.72 * Pr^0.37 * (t.pitch / d_hyd)^(-0.372) : 3.66
    return Correlations.get_h_from_Nu(Nu, k, d_hyd)
end

"""Custom (user-defined power-law) turbulator."""
mutable struct CustomTurbulator
    turbulator_type::String
    a::Union{Float64,Nothing}
    b::Union{Float64,Nothing}
    c::Union{Float64,Nothing}
end
CustomTurbulator(; a=nothing, b=nothing, c=nothing) =
    CustomTurbulator("Custom", a, b, c)

function k_t_correlation(t::CustomTurbulator;
                         Re::Float64, Sc::Float64,
                         d_hyd::Float64, D::Float64)::Float64
    Sh = Re > 2030 ? t.a * Re^t.b * Sc^t.c : 3.66
    return Correlations.get_k_from_Sh(Sh, d_hyd, D)
end

function h_t_correlation(t::CustomTurbulator;
                         Re::Float64, Pr::Float64,
                         d_hyd::Float64, k::Float64)::Float64
    Nu = Re > 2030 ? t.a * Re^t.b * Pr^t.c : 3.66
    # Note: the Python source references self.pitch here, which appears to be a
    # copy-paste bug (CustomTurbulator has no pitch). Using d_hyd as characteristic
    # length, consistent with the WireCoil implementation.
    return Correlations.get_h_from_Nu(Nu, k, d_hyd)
end

# ── Cross-sections ─────────────────────────────────────────────────────────────
# The flow cross-section of a pipe. Each concrete type knows its flow area and
# wetted perimeter, from which the hydraulic diameter D_h = 4 A / P follows.
# Circular is the default; non-circular shapes (Rectangular, Annulus, …) plug in
# here without touching the transport code, which only ever sees `d_Hyd`.

"""Abstract supertype for pipe flow cross-sections."""
abstract type CrossSection end

"""Circular cross-section of inner diameter `D` [m] (the default pipe shape)."""
struct Circular <: CrossSection
    D::Float64
end

"""Rectangular cross-section of `width` × `height` [m]."""
struct Rectangular <: CrossSection
    width::Float64
    height::Float64
end

"""Concentric annular cross-section between `D_outer` and `D_inner` [m]."""
struct Annulus <: CrossSection
    D_outer::Float64
    D_inner::Float64
end

"""
    TwistedElliptical(a, b)

Elliptical cross-section with semi-major axis `a` and semi-minor axis `b` [m].
Intended for twisted-elliptical tubes (TETs): the twist pitch affects swirl-induced
heat/mass-transfer enhancement but not the geometric `D_h = 4A/P`, which depends
only on the cross-sectional area and perimeter.
"""
struct TwistedElliptical <: CrossSection
    a::Float64   # semi-major axis [m]
    b::Float64   # semi-minor axis [m]
end

"""Flow (wetted) cross-sectional area [m²]."""
flow_area(cs::Circular)::Float64         = π * (cs.D / 2)^2
flow_area(cs::Rectangular)::Float64      = cs.width * cs.height
flow_area(cs::Annulus)::Float64          = π * ((cs.D_outer / 2)^2 - (cs.D_inner / 2)^2)
flow_area(cs::TwistedElliptical)::Float64 = π * cs.a * cs.b

"""Wetted perimeter [m]."""
wetted_perimeter(cs::Circular)::Float64    = π * cs.D
wetted_perimeter(cs::Rectangular)::Float64 = 2 * (cs.width + cs.height)
wetted_perimeter(cs::Annulus)::Float64     = π * (cs.D_outer + cs.D_inner)

# Ramanujan's second approximation; exact for circles (h=0), O(h^10) error otherwise.
function wetted_perimeter(cs::TwistedElliptical)::Float64
    a, b = cs.a, cs.b
    h = ((a - b) / (a + b))^2
    return π * (a + b) * (1 + 3h / (10 + sqrt(4 - 3h)))
end

"""Hydraulic diameter of a cross-section, ``D_h = 4 A / P``."""
hydraulic_diameter(cs::CrossSection)::Float64 =
    hydraulic_diameter(flow_area(cs), wetted_perimeter(cs))

# ── Geometry ───────────────────────────────────────────────────────────────────

"""
    Geometry

Geometric parameters for a pipe/tube component.

Fields
- `L`           : length [m]
- `D`           : inner diameter [m] (circular pipes; characteristic diameter otherwise)
- `ds`          : wall thickness [m]
- `n_pipes`     : number of parallel pipes (default 1)
- `turbulator`  : optional turbulator object
- `cross_section`: optional `CrossSection`. When `nothing` (default) the pipe is
  treated as `Circular(D)`, so the hydraulic diameter equals `D`. Provide a
  `Rectangular`, `Annulus`, … to model non-cylindrical channels.
"""
mutable struct Geometry
    L::Union{Float64,Nothing}
    D::Union{Float64,Nothing}
    ds::Union{Float64,Nothing}
    n_pipes::Float64
    turbulator::Union{WireCoil,CustomTurbulator,Turbulator,Nothing}
    cross_section::Union{CrossSection,Nothing}
end
function Geometry(; L=nothing, D=nothing, ds=nothing, n_pipes=1.0,
                    turbulator=nothing, cross_section=nothing)
    Geometry(L, D, ds, n_pipes, turbulator, cross_section)
end

# Effective cross-section: the explicit one if given, else a circular pipe of
# diameter `D`. Built on demand so it always reflects the current `D`.
_section(g::Geometry)::CrossSection =
    g.cross_section !== nothing ? g.cross_section : Circular(g.D)

"""Flow (wetted) cross-sectional area of the pipe [m²]."""
flow_area(g::Geometry)::Float64 = flow_area(_section(g))

"""Wetted perimeter of the pipe cross-section [m]."""
wetted_perimeter(g::Geometry)::Float64 = wetted_perimeter(_section(g))

"""Hydraulic diameter of the pipe, ``D_h = 4 A / P`` (equals `D` for a circular pipe)."""
hydraulic_diameter(g::Geometry)::Float64 = hydraulic_diameter(_section(g))

"""Fluid volume of a single pipe [m³]."""
get_fluid_volume(g::Geometry)::Float64 = flow_area(g) * g.L

"""
Wall (solid) volume of a single pipe [m³].

Currently modeled as a circular wall of thickness `ds` (inner diameter `D`),
so a circular `D` is required even when a non-circular `cross_section` is set.
"""
function get_solid_volume(g::Geometry)::Float64
    g.D === nothing && error(
        "get_solid_volume requires a circular diameter `D`; wall volume for " *
        "non-circular cross-sections is not yet modeled")
    return π * ((g.D / 2)^2 - (g.D / 2 - g.ds)^2) * g.L
end

"""Total volume (fluid + wall) of a single pipe [m³]."""
get_total_volume(g::Geometry)::Float64 = get_fluid_volume(g) + get_solid_volume(g)

# ── FluidMaterial / SolidMaterial (simple property containers) ─────────────────

"""Physical properties of a fluid material at a given temperature."""
mutable struct FluidMaterial
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
mutable struct SolidMaterial
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
mutable struct Fluid
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
function set_properties_from_fluid_material!(fluid::Fluid, mat::FluidMaterial)::Nothing
    fluid.T         = mat.T
    fluid.D         = mat.D
    fluid.Solubility = mat.Solubility
    fluid.mu        = mat.mu
    fluid.rho       = mat.rho
    fluid.cp        = mat.cp
    fluid.k         = mat.k
    return nothing
end

"""Re-evaluate temperature-dependent D and Solubility using Arrhenius expressions."""
function update_T_prop!(fluid::Fluid)::Nothing
    if fluid.D_0 !== nothing && fluid.E_d !== nothing && fluid.T !== nothing
        fluid.D = fluid.D_0 * exp(-fluid.E_d / (K_B_EV * fluid.T))
    end
    if fluid.Solubility_0 !== nothing && fluid.E_s !== nothing && fluid.T !== nothing
        fluid.Solubility = fluid.Solubility_0 * exp(-fluid.E_s / (K_B_EV * fluid.T))
    end
    return nothing
end

"""
    get_kt!(fluid::Fluid; turbulator=nothing)

Compute the mass-transfer coefficient `k_t` from the fluid's hydraulic diameter,
velocity, and thermophysical properties. Uses the Getthem correlation for turbulent
flow and Sh = 3.66 for laminar flow, or the turbulator's own correlation if one is
provided.
"""
function get_kt!(fluid::Fluid;
                 turbulator::Union{WireCoil,CustomTurbulator,Turbulator,Nothing}=nothing)::Nothing
    if fluid.d_Hyd === nothing
        @warn "Hydraulic diameter is not defined"
        return nothing
    end
    if fluid.k_t !== nothing
        println("k_t is already defined")
        return nothing
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
    return nothing
end

# ── Membrane ───────────────────────────────────────────────────────────────────

"""
    Membrane

Metallic membrane properties for hydrogen-isotope transport.

Temperature-dependent `D` and `K_S` are evaluated from Arrhenius parameters
(`D_0`/`E_d` and `K_S_0`/`E_S`) if provided.
"""
mutable struct Membrane
    T::Union{Float64,Nothing}
    D::Union{Float64,Nothing}
    ds::Union{Float64,Nothing}
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
    T=nothing, D=nothing, ds=nothing, K_S=nothing,
    k_d=nothing, k_r=nothing, k=nothing,
    D_0=nothing, E_d=nothing, K_S_0=nothing, E_S=nothing,
    inv=nothing, V=nothing
)
    D_val  = (D_0 !== nothing && E_d !== nothing && T !== nothing) ?
        D_0 * exp(-E_d / (K_B_EV * T)) : D
    KS_val = (K_S_0 !== nothing && E_S !== nothing && T !== nothing) ?
        K_S_0 * exp(-E_S / (K_B_EV * T)) : K_S
    Membrane(T, D_val, ds, KS_val, k_d, k_r, k, D_0, E_d, K_S_0, E_S, inv, V)
end

"""Copy membrane properties from a `SolidMaterial` object."""
function set_properties_from_solid_material!(mem::Membrane, mat::SolidMaterial)::Nothing
    mem.T   = mat.T
    mem.D   = mat.D
    mem.K_S = mat.K_S
    return nothing
end

"""Re-evaluate temperature-dependent D and K_S using Arrhenius expressions."""
function update_T_prop!(mem::Membrane)::Nothing
    if mem.D_0 !== nothing && mem.E_d !== nothing && mem.T !== nothing
        mem.D = mem.D_0 * exp(-mem.E_d / (K_B_EV * mem.T))
    end
    if mem.K_S_0 !== nothing && mem.E_S !== nothing && mem.T !== nothing
        mem.K_S = mem.K_S_0 * exp(-mem.E_S / (K_B_EV * mem.T))
    end
    return nothing
end

end # module PipeSubclasses
