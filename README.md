# TRIOMA.jl

**TRItium Object-oriented and Modular Analysis — Julia port**

A Julia implementation of the [TRIOMA Python package](https://github.com/gabriele-ferrero/TRIOMA) by Gabriele Ferrero and Samuele Meschini (Politecnico di Torino). TRIOMA provides physics-based models for tritium transport in fusion reactor outer fuel cycles (OFC): extraction efficiency, inventories, and losses in individual components and full circuits.

---

## Features

- **Component types:** Permeation Against Vacuum (PAV) pipes, Gas-Liquid Contactors (packed columns), Breeding Blankets, Heat Exchangers
- **Fluid types:** Molten salts (Henry's law) and liquid metals (Sievert's law), with Arrhenius-evaluated diffusivity and solubility
- **Pre-built materials:** Flibe, Sodium, LiPb fluids; Steel membranes
- **Transport regimes:** Automatic identification of mass-transport-limited, diffusion-limited, surface-limited, and mixed regimes
- **Efficiency methods:** Analytical (Lambert W closed-form) and numerical (1-D integration) extraction efficiency
- **Inventory:** Fluid-phase and membrane-phase tritium inventory, analytical and numerical
- **Circuit analysis:** Connect components in series, solve closed loops to steady state, compute net gains/losses and pumping power
- **Turbulators:** Wire-coil and user-defined power-law turbulator correlations
- **Non-circular cross-sections:** `Circular`, `Rectangular`, `Annulus`, and `TwistedElliptical` cross-section types with hydraulic diameter `D_h = 4A/P` computed automatically
- **Twisted-tube correlations:** Tube-side Nusselt and friction-factor correlations for twisted elliptical tubes (Ievlev, Asmantas, Simonin, Yang, Gao) verified against Hughes (2017)

---

## Installation

TRIOMA.jl is not yet registered in the General registry. Install from a local path:

```julia
using Pkg
Pkg.add(path="/path/to/TRIOMA_JL")
```

**Runtime dependencies** (resolved automatically by Pkg):

| Package | Purpose |
|---|---|
| `AtomicAndPhysicalConstants` | Boltzmann constant for Arrhenius evaluations |
| `LambertW` | Lambert W function for the analytical molten-salt efficiency formula |
| `Optim` | Nonlinear minimization for mixed-regime flux and c_out solvers |
| `QuadGK` | Adaptive Gauss-Kronrod quadrature for NTU integrals and inventory |

Julia **≥ 1.9** is required.

---

## Package structure

```
TRIOMA_JL/
├── Project.toml
├── README.md
├── src/
│   ├── TRIOMA.jl               # Top-level module; loads submodules, re-exports public API
│   ├── TriomaCore.jl           # TriomaClass struct, inspect, update_attribute!
│   ├── Correlations.jl         # Nu, Sh, Re, Pr, f correlations; twisted-tube correlations
│   ├── PipeSubclasses.jl       # Geometry, Fluid, Membrane, FluidMaterial, SolidMaterial,
│   │                           #   CrossSection hierarchy, turbulator types
│   ├── Materials.jl            # Flibe, Sodium, LiPb, Steel property functions
│   ├── RegimeHelpers.jl        # FusionCoolant module: H/W dimensionless numbers,
│   │                           #   regime classifiers for molten salts and liquid metals
│   ├── ExtractorFunctions.jl   # NTU integrals and c_out solvers for packed columns
│   ├── GasLiquidContactor.jl   # GLC_Gas and GLC structs + methods
│   ├── PAV.jl                  # Component struct with all efficiency/flux/inventory methods
│   ├── BreedingBlanket.jl      # BreedingBlanket struct and get_cout!
│   └── Circuit.jl              # Circuit struct, solve_circuit!, get_eff_circuit!, etc.
└── test/
    ├── runtests.jl
    ├── test_correlations.jl
    ├── test_trioma_class.jl
    ├── test_component_tools.jl
    ├── test_hydraulic_diameter.jl
    ├── test_imports.jl
    └── test_project.jl
```

A single `using TRIOMA` brings the entire public API into scope — no submodule imports needed.

---

## Quick start

### 1. Single PAV pipe — Flibe molten salt, analytical efficiency

```julia
using TRIOMA

T     = 973.15   # K — operating temperature
d_hyd = 25.4e-3  # m — inner diameter (circular pipe)
U0    = 2.5      # m/s — bulk fluid velocity

# Build fluid from Flibe material properties at T
flibe_mat = Flibe(T)
flibe = Fluid(d_Hyd=d_hyd, U0=U0, MS=true)
set_properties_from_fluid_material!(flibe, flibe_mat)

# Steel membrane
steel_mat = Steel(T)
steel = Membrane(dw=0.25e-3, k_r=1e9, k_d=1e9)
set_properties_from_solid_material!(steel, steel_mat)

# Geometry and component
geom = Geometry(L=1.0, D=d_hyd, dw=0.25e-3)
pav  = Component(c_in=1e-3, geometry=geom, fluid=flibe, membrane=steel, name="PAV")

# Volumes and hydraulic diameter (auto-derives d_Hyd = D for circular pipe)
define_component_volumes!(pav)

# Compute mass-transfer coefficient and extraction efficiency
get_kt!(pav.fluid)
use_analytical_efficiency!(pav)   # Lambert W closed-form
outlet_c_comp!(pav)

println("Efficiency : ", round(pav.eff * 100, digits=2), " %")
println("c_out      : ", pav.c_out, " mol/m³")
```

### 2. Transport regime identification

```julia
regime = get_regime(pav; print_var=true)
# Prints H, W, H/W and the regime label, e.g.:
#   H = 0.0234
#   W = 1.2e-5
#   H/W = 1950.0
#   Mass transport limited
```

### 3. Inspect and update fields

```julia
inspect(pav)                          # print all fields recursively
inspect(pav; variable_name="eff")    # print one field by name

update_attribute!(pav, "c_in", 5e-3) # update any field by name string
use_analytical_efficiency!(pav)
outlet_c_comp!(pav)
```

### 4. Liquid-metal component (Sodium, Sievert's law)

```julia
T_na   = 750.0
sodium = Fluid(d_Hyd=20e-3, U0=1.5, MS=false)
set_properties_from_fluid_material!(sodium, Sodium(T_na))

steel_na = Membrane(dw=3e-4, k_r=1e8, k_d=1e8)
set_properties_from_solid_material!(steel_na, Steel(T_na))

geom_na = Geometry(L=2.0, D=20e-3, dw=3e-4)
comp_na = Component(c_in=2e-4, geometry=geom_na, fluid=sodium, membrane=steel_na)
define_component_volumes!(comp_na)
get_kt!(comp_na.fluid)
use_analytical_efficiency!(comp_na)
println("LM efficiency: ", round(comp_na.eff * 100, digits=2), " %")
```

### 5. Non-circular tube — TwistedElliptical cross-section

```julia
# Semi-major axis a = 11 mm, semi-minor axis b = 4.5 mm, twist pitch s = 0.2 m
cs   = TwistedElliptical(0.011, 0.0045)
geom = Geometry(L=1.0, D=0.022, dw=5e-4, cross_section=cs)

# Hydraulic diameter is derived from 4A/P automatically
define_component_volumes!(comp)    # sets comp.fluid.d_Hyd = hydraulic_diameter(cs)

println("D_h = ", hydraulic_diameter(cs), " m")
println("Flow area = ", flow_area(cs), " m²")

# Tube-side Nusselt number (Simonin correlation, valid 1000 ≤ Re ≤ 17 000)
D_max = 2 * cs.a
D_min = 2 * cs.b
D_h   = hydraulic_diameter(cs)
s     = 0.2  # twist pitch [m]
Re_t  = 8000.0
Pr_t  = 12.0
Nu    = Nu_Si_tube(Re_t, Pr_t, D_h, D_max, s)
println("Nu (Simonin) = ", round(Nu, digits=1))
```

### 6. Gas-Liquid Contactor (packed column)

```julia
T_glc = 673.0          # K
K_S   = 1.33e-4 * exp(-1350.0 / (8.314 * T_glc))  # Sievert constant

sweep = GLC_Gas(G_gas=3e-3/3600, pg_in=0.0, p_tot=1.5e5)
pbli  = Fluid(Solubility=K_S, MS=false)

glc = GLC(
    H=0.6, R=0.0547/2,
    c_in=1085.0^0.5 * K_S, c_out=653.0^0.5 * K_S,
    fluid=pbli, GLC_gas=sweep, T=T_glc, G_L=71e-3/3600,
)

B_l, kla = get_kla_from_cout!(glc)
println("kla = ", round(kla * 1e3, digits=3), " × 10⁻³ s⁻¹")
```

### 7. Circuit

```julia
# Connect components and solve a closed loop to steady state
circ = Circuit()
add_component!(circ, pav)
add_component!(circ, comp_na)
solve_circuit!(circ)
get_eff_circuit!(circ)
inspect_circuit(circ)
```

---

## API reference

### Types

| Type | Description |
|---|---|
| `Geometry` | Pipe dimensions: `L`, `D`, `dw` (wall thickness), `n_pipes`, `turbulator`, `cross_section` |
| `Fluid` | Fluid transport properties; Arrhenius evaluation of `D`/`Solubility` if `D_0`/`E_d` given |
| `Membrane` | Metal wall properties; Arrhenius evaluation of `D`, `K_S` |
| `FluidMaterial` | Property container returned by `Flibe(T)`, `Sodium(T)`, `LiPb(T)` |
| `SolidMaterial` | Property container returned by `Steel(T)` |
| `WireCoil` | Wire-coil turbulator; pitch-dependent `k_t` and `h` correlations |
| `CustomTurbulator` | User-defined Sh = a·Reᵇ·Scᶜ turbulator |
| `Component` | Full PAV/HX permeation component |
| `GLC_Gas` | Sweep-gas side of a packed column |
| `GLC` | Gas-Liquid Contactor packed column |
| `BreedingBlanket` | Tritium-producing breeding blanket |
| `Circuit` | Series circuit of `Component` and `BreedingBlanket` objects |
| `TriomaClass` | Generic inspectable/updatable container (testing and scripting) |

### Cross-section types

All live in the `CrossSection` abstract hierarchy. Constructors take inner dimensions.

| Type | Constructor | `D_h = 4A/P` |
|---|---|---|
| `Circular` | `Circular(D)` | equals `D` |
| `Rectangular` | `Rectangular(width, height)` | `2wh/(w+h)` |
| `Annulus` | `Annulus(D_outer, D_inner)` | `D_outer − D_inner` |
| `TwistedElliptical` | `TwistedElliptical(a, b)` | `4πab / P` (Ramanujan perimeter) |

`flow_area`, `wetted_perimeter`, and `hydraulic_diameter` dispatch on any `CrossSection`.
Set `Geometry(cross_section=cs)` to override the default circular assumption;
`define_component_volumes!` then derives `d_Hyd` automatically.

### Core functions — Component

| Function | Description |
|---|---|
| `define_component_volumes!(comp)` | Set fluid/membrane volumes and auto-derive `d_Hyd` from geometry |
| `get_kt!(fluid)` | Compute mass-transfer coefficient `k_t` (Getthem or turbulator correlation) |
| `set_hydraulic_diameter!(comp)` | Explicitly derive and store `d_Hyd = 4A/P` from geometry |
| `use_analytical_efficiency!(comp)` | Compute efficiency via Lambert W / closed-form; store in `comp.eff` |
| `analytical_efficiency!(comp)` | As above but fills `comp.eff_an` only |
| `get_efficiency!(comp)` | Numerical 1-D integration along pipe |
| `outlet_c_comp!(comp)` | Fill `comp.c_out`, accounting for recirculation |
| `get_regime(comp)` | Return regime string; optionally print H, W, H/W |
| `get_adimensionals!(comp)` | Compute and store H and W |
| `get_inventory!(comp)` | Total tritium inventory (fluid + membrane) |
| `get_solid_inventory!(comp)` | Membrane inventory only |
| `get_fluid_inventory!(comp)` | Fluid inventory only |
| `get_pressure_drop!(comp)` | Darcy-Weisbach pressure drop |
| `get_pumping_power!(comp)` | Pumping power [W] |
| `get_global_HX_coeff!(comp)` | Overall heat transfer coefficient U |
| `inspect(obj)` | Print all or named fields recursively |
| `update_attribute!(obj, name, val)` | Recursively set a field by name string |

### Material functions

| Call | Returns | Notes |
|---|---|---|
| `Flibe(T)` | `FluidMaterial` | Molten salt; `MS=true` |
| `Sodium(T)` | `FluidMaterial` | Liquid metal; `MS=false` |
| `LiPb(T)` | `FluidMaterial` | Liquid metal; some properties are TODO placeholders |
| `Steel(T)` | `SolidMaterial` | `K_S = 1.0` is a placeholder |

### Twisted-tube (TET) correlations

All take scalar `Float64` arguments. `D_max = 2a`, `D_min = 2b`, `D_h = hydraulic_diameter(cs)`,
`s` = twist pitch in metres (axial length per full 360° rotation).

**Nusselt number — tube-side:**

| Function | Correlation | Valid Re range | Notes |
|---|---|---|---|
| `Nu_Iev_tube(Re, D_max, s)` | Ievlev | 6 000 – 100 000 | No Pr dependence |
| `Nu_Asma_tube(Re, Pr, D_max, s)` | Asmantas Modified | 7 000 – 200 000 | ϕ = 1 (no wall correction) |
| `Nu_Si_tube(Re, Pr, D_h, D_max, s)` | Simonin | 1 000 – 17 000 | — |
| `Nu_Yang_lam(Re, Pr, D_max, D_min, s)` | Yang Laminar | 100 – 500 | Best overall fit per Hughes (2017) |
| `Nu_Yang1_tube(Re, Pr, D_h, D_max, D_min, s)` | Yang #1 | 5 000 – 20 000 | Narrow geometry range |
| `Nu_Yang2_tube(Re, Pr, D_h, D_max, D_min, s)` | Yang #2 | 7 900 – 26 500 | Pr in denominator (as published) |

**Darcy friction factor — tube-side:**

| Function | Correlation | Valid Re range |
|---|---|---|
| `fD_Iev_tube(Re, D_max, s)` | Ievlev | 6 000 – 100 000 |
| `fD_Asma_tube(Re, D_max, s)` | Asmantas | 7 000 – 200 000 |
| `fD_Si_tube(Re, D_max, s)` | Simonin | 1 000 – 17 000 |
| `fD_Gao_tube(Re, D_h, D_max, D_min, s)` | Gao | 5 000 – 20 000 |
| `fD_Yang_tube(Re, Pr, D_h, D_max, D_min, s)` | Yang | 7 900 – 26 500 |

All correlations verified against Appendix A of Hughes, J.T. (2017), *Experimental and Computational Investigations of Heat Transfer Systems in Fluoride Salt-cooled High-temperature Reactors*, University of New Mexico PhD dissertation.

### General correlations (`Correlations` module)

| Function | Description |
|---|---|
| `Nu_DittusBoelter(Re, Pr)` | Standard turbulent pipe flow Nusselt number |
| `Nu_Gnielinsky(Re, Pr, f)` | Gnielinski correlation |
| `Nu_SiederTate(Re, Pr, mu_w, mu_c)` | Sieder-Tate with viscosity correction |
| `Re(rho, u, L, mu)` | Reynolds number |
| `Pr(c_p, mu, k)` | Prandtl number |
| `Schmidt(D, mu, rho)` | Schmidt number |
| `Sherwood(Sc, Re)` | Getthem correlation (default in `get_kt!`) |
| `Sherwood_HT_analogy(Re, Sc)` | Heat-transfer analogy Sherwood number |
| `hydraulic_diameter(area, perimeter)` | `D_h = 4A/P` (scalar primitive) |
| `get_h_from_Nu(Nu, k, D)` | Heat transfer coefficient from Nu |
| `get_k_from_Sh(Sh, L, D)` | Mass transfer coefficient from Sh |
| `f_Haaland(Re, e_D)` | Haaland friction factor |
| `get_length_HX(dTML, d_hyd, U, Q)` | HX tube length from LMTD |
| `get_deltaTML(...)` | Log-mean temperature difference |

---

## Differences from the Python version

| Topic | Python | Julia |
|---|---|---|
| Imports | `from TRIOMA import Component, ...` | `using TRIOMA` (all exports available) |
| Mutation style | Methods mutate `self` | Functions with `!` suffix mutate first argument |
| Wall thickness field | `thick` | `dw` (in `Geometry` and `Membrane`) |
| Module layout | Nested package with `tools/` subdirectory | Flat `src/` with `include`-based submodules |
| Regime helpers | `MoltenSalts`, `LiquidMetals` modules | Unified `FusionCoolant` module |
| Plotting | `plot_component()`, `plot_circuit()` via matplotlib | Not ported |
| `scipy.integrate.nquad` | 2-D inventory integrals | Nested `quadgk` calls |
| `scipy.optimize.minimize` (Powell) | Flux/c_out solver | `Optim.optimize` with `NelderMead()` |
| Lambert W | `scipy.special.lambertw` | `LambertW.jl` |

---

## Known limitations

- **Plotting** (`plot_component`, `plot_circuit`) is not ported.
- **`LiPb`** thermal conductivity and specific heat are placeholder values (`1.0`) — marked `# TODO` in the upstream Python code as well.
- **`Steel`** Sievert constant `K_S = 1.0` is a placeholder in both versions.
- **`get_solid_volume`** uses the circular-wall formula `π((D/2)² − (D/2 − dw)²)L` regardless of `cross_section`. Wall volume for non-circular cross-sections is not yet modelled.

---

## Citations

Correlations and physical models are drawn from:

- Alberghi, C. et al. "Development of new analytical tools for tritium transport modelling." *Fusion Engineering and Design* 177 (2022): 113083.
- Humrickhouse, P. W. and T. F. Fuerst. *Tritium transport phenomena in molten-salt reactors.* INL/EXT-20-59927, Idaho National Lab., 2020.
- Fuerst, T. F., C. N. Taylor and P. W. Humrickhouse. *Molten Salt Tritium Transport Experiment Design.* INL/EXT-21-63108, Idaho National Lab., 2021.
- Rader, J. D. et al. "Verification of modelica-based models with analytical solutions for tritium diffusion." *Nuclear Technology* 203.1 (2018): 58–65.
- Hughes, J. T. *Experimental and Computational Investigations of Heat Transfer Systems in Fluoride Salt-cooled High-temperature Reactors.* PhD dissertation, University of New Mexico, 2017. *(Source for twisted-tube TET correlations.)*
- Dzyubenko, B. V., Dreitser, G. A., Ashmantas, L-V. A. and Segal, M. D. *Impaired and Enhanced Heat Transfer in Helically Finned Tubes and Tube Bundles.* Hemisphere, Washington, 1990.

## License

MIT — same as the upstream TRIOMA Python package.
