# TRIOMA.jl

**TRItium Object-oriented and Modular Analysis — Julia port**

A Julia translation of the [TRIOMA Python package](https://github.com/gabriele-ferrero/TRIOMA) by Gabriele Ferrero and Samuele Meschini (Politecnico di Torino). TRIOMA helps engineers design the Outer Fuel Cycle (OFC) of fusion reactors by providing compact, physics-based functions for estimating tritium extraction efficiency, losses, and inventories in individual components and full circuits.

> **Accuracy note:** This translation follows the original Python source closely, including equations, correlations, and numerical methods. A small number of material properties (LiPb properties, Steel Sievert constant, Sodium thermal properties) carry `TODO` placeholders that are also present in the upstream Python code — these are noted inline.

---

## Features

- **Component types:** Permeation Against Vacuum (PAV) pipes, Gas-Liquid Contactors (packed columns), Breeding Blankets, Heat Exchangers
- **Fluid types:** Molten salts (Henry's law) and liquid metals (Sievert's law)
- **Pre-built materials:** Flibe, Sodium, LiPb fluids; Steel membranes
- **Transport regimes:** Automatic identification of mass-transport-limited, diffusion-limited, surface-limited, and mixed regimes
- **Efficiency methods:** Analytical (Lambert W / closed-form) and numerical (1-D integration) extraction efficiency
- **Inventory:** Fluid-phase and membrane-phase tritium inventory, analytical and numerical
- **Circuit analysis:** Connect components in series, solve closed loops to steady state, compute gains/losses/pumping power
- **Turbulators:** Wire coil and user-defined custom turbulator correlations

---

## Installation

This package is not yet registered in the Julia General registry. Install directly from the local path (or a future GitHub URL):

```julia
# From Julia's package manager (Pkg REPL — press ] to enter)
pkg> add /path/to/TRIOMA_jl

# Or from the Julia REPL
using Pkg
Pkg.add(path="/path/to/TRIOMA_jl")
```

**Dependencies** (resolved automatically by Pkg):

| Package | Purpose |
|---|---|
| `Optim` | Nonlinear minimization for mixed-regime flux and c_out solvers |
| `QuadGK` | Adaptive Gauss-Kronrod quadrature for NTU integrals and inventory |
| `LambertW` | Lambert W function for the analytical molten-salt efficiency formula |
| `SpecialFunctions` | Supporting special functions |

Julia **≥ 1.9** is required (uses `@match`-style dispatch patterns via `if/elseif` chains).

---

## Package structure

```
TRIOMA_jl/
├── Project.toml
├── README.md
└── src/
    └── tools/
        ├── TriomaClass.jl          # Abstract base type, inspect, update_attribute!
        ├── Correlations.jl         # Nu, Sh, Re, Pr, friction-factor correlations
        ├── Materials.jl            # Flibe, Sodium, LiPb, Steel property functions
        ├── MoltenSalts.jl          # H_ms, W_ms dimensionless numbers, regime classifier
        ├── LiquidMetals.jl         # W_lm, partition_param, regime classifier
        └── Extractors/
            ├── PipeSubclasses.jl   # Geometry, Fluid, Membrane, FluidMaterial,
            │                       #   SolidMaterial, WireCoil, CustomTurbulator
            ├── Extractor.jl        # NTU integrals and c_out solvers for packed columns
            ├── GasLiquidContactor.jl  # GLC_Gas and GLC structs + methods
            ├── PAV.jl              # Component struct with all efficiency/flux/inventory methods
            ├── BreedingBlanket.jl  # BreedingBlanket struct and get_cout!
            └── Circuit.jl          # Circuit struct, solve_circuit!, get_eff_circuit!, etc.
```

> **Note:** `BreedingBlanket.jl` and `Circuit.jl` are listed here as part of the complete package specification but may need to be added if not already present — see the Python originals at `src/TRIOMA/tools/BreedingBlanket.py` and `Circuit.py` for the full reference implementation.

---

## Quick start

### 1. Single PAV pipe — analytical efficiency (Flibe / molten salt)

This reproduces the core example from the original Python tutorial.

```julia
using TRIOMA.TriomaTypes
using TRIOMA.PipeSubclasses
using TRIOMA.Materials
using TRIOMA.PAVModule

# --- Operating conditions ---
T    = 973.15   # K
d_hyd = 25.4e-3  # m  (inner diameter)
U0    = 2.5      # m/s

# --- Fluid: Flibe molten salt ---
flibe_mat = Flibe(T)                        # FluidMaterial at temperature T
flibe     = Fluid(d_Hyd=d_hyd, U0=U0, MS=true)
set_properties_from_fluid_material!(flibe, flibe_mat)

# --- Membrane: Steel wall ---
steel_mat = Steel(T)
steel     = Membrane(thick=0.25e-3, k_r=1e9, k_d=1e9)
set_properties_from_solid_material!(steel, steel_mat)

# --- Geometry and component ---
geom  = Geometry(L=1.0, D=25.4e-3, thick=0.25e-3)
c_in  = 1e-3                                # mol/m³ inlet concentration
pav   = Component(c_in=c_in, geometry=geom, fluid=flibe, membrane=steel, name="PAV")

# --- Compute efficiency ---
use_analytical_efficiency!(pav)             # fills pav.eff via Lambert W formula
outlet_c_comp!(pav)                         # fills pav.c_out = c_in * (1 - eff)

println("Extraction efficiency : ", round(pav.eff * 100, digits=2), " %")
println("Outlet concentration  : ", pav.c_out, " mol/m³")
```

### 2. Inspect a component

```julia
# Print all fields (recurses into nested Fluid and Membrane)
inspect(pav)

# Print only one field by name (case-insensitive)
inspect(pav; variable_name="eff")
inspect(pav; variable_name="eff_an")

println("Analytical vs numerical relative error: ",
        abs(pav.eff - pav.eff_an) / pav.eff)
```

### 3. Update a field and recompute

```julia
# Change inlet concentration and rerun
update_attribute!(pav, "c_in", 5e-3)
use_analytical_efficiency!(pav)
outlet_c_comp!(pav)
println("New efficiency at c_in=5e-3: ", round(pav.eff * 100, digits=2), " %")
```

### 4. Identify the transport regime

```julia
regime = get_regime(pav; print_var=true)
# Prints H, W, H/W and the regime label, e.g.:
#   H = 0.0234
#   W = 1.2e-5
#   H/W = 1950.0
#   Mass transport limited
```

### 5. Tritium inventory

```julia
# Analytical membrane + fluid inventory
get_inventory!(pav; flag_an=true)
println("Membrane inventory : ", pav.membrane.inv, " mol")
println("Fluid inventory    : ", pav.fluid.inv,    " mol")
println("Total inventory    : ", pav.inv,           " mol")
```

### 6. Numerical efficiency (full 1-D integration along pipe)

```julia
# Slower but does not require the analytical regime assumptions
get_efficiency!(pav; c_guess=c_in * 1e-3)
println("Numerical efficiency: ", round(pav.eff * 100, digits=2), " %")
```

### 7. Liquid-metal component (Sodium / Sievert's law)

```julia
T_na   = 750.0
Na_mat = Sodium(T_na)
sodium = Fluid(d_Hyd=20e-3, U0=1.5, MS=false)
set_properties_from_fluid_material!(sodium, Na_mat)

steel_na = Membrane(thick=3e-4, k_r=1e8, k_d=1e8)
set_properties_from_solid_material!(steel_na, Steel(T_na))

geom_na = Geometry(L=2.0, D=20e-3, thick=3e-4)
comp_na = Component(c_in=2e-4, geometry=geom_na, fluid=sodium, membrane=steel_na, name="Na_pipe")

use_analytical_efficiency!(comp_na)
println("LM pipe efficiency: ", round(comp_na.eff * 100, digits=2), " %")
```

### 8. Gas-Liquid Contactor (packed column)

```julia
using TRIOMA.GasLiquidContactorModule

R_const = 8.314
T_glc   = 673.0          # K
Z       = 0.6            # m column height
R_col   = 0.0547 / 2     # m column radius
G_l     = 71e-3 / 3600   # m³/s liquid flow rate
G_gas   = 3e-3  / 3600   # m³/s gas flow rate
p_in    = 1085.0         # Pa tritium partial pressure inlet
p_out_l = 653.0          # Pa tritium partial pressure outlet (known, for kla back-calc)

# PbLi solubility (Sievert)
K_S = 1.33e-4 * exp(-1350.0 / (R_const * T_glc))

sweep_gas = GLC_Gas(G_gas=G_gas, pg_in=0.0, p_tot=1.5e5)
pb_li     = Fluid(Solubility=K_S, MS=false)

melodie = GLC(
    H       = Z,
    R       = R_col,
    c_in    = p_in^0.5  * K_S,
    c_out   = p_out_l^0.5 * K_S,  # known outlet for back-calculation
    fluid   = pb_li,
    GLC_gas = sweep_gas,
    T       = T_glc,
    G_L     = G_l,
)

B_l, kla = get_kla_from_cout!(melodie)
println("Liquid holdup B_l : ", round(B_l, digits=2), " m/h")
println("kla               : ", round(kla * 1e3, digits=3), " × 10⁻³ s⁻¹")

# Forward: given kla, compute c_out and efficiency
update_attribute!(melodie, "c_out", nothing)   # reset
melodie.kla = kla
get_c_out!(melodie)
println("Computed efficiency: ", round(melodie.eff * 100, digits=2), " %")
```

### 9. Wire-coil turbulator

```julia
wc   = WireCoil(pitch=5e-3)
geom = Geometry(L=1.5, D=25.4e-3, thick=0.25e-3, turbulator=wc)
# Fluid k_t will be computed using the WireCoil correlation when needed
pav_wc = Component(c_in=1e-3, geometry=geom, fluid=flibe, membrane=steel)
use_analytical_efficiency!(pav_wc)
println("Wire-coil efficiency: ", round(pav_wc.eff * 100, digits=2), " %")
```

### 10. Heat exchanger — temperature-discretized circuit

```julia
# A Component with heat-exchanger global U-value can be split into N
# sub-components that capture the temperature profile along the tube.
# Worth verifying this against the HX tutorial notebook in the Python source
# for your specific operating conditions.

get_global_HX_coeff!(pav)   # fills pav.U
# Then use pav.split_HX(...) — see PAV.jl for the full signature
```

---

## Key API reference

### Types

| Julia type | Python equivalent | Description |
|---|---|---|
| `Geometry` | `Geometry` | Pipe dimensions (L, D, thick, n_pipes, turbulator) |
| `Fluid` | `Fluid` | Fluid transport properties; computes D/Solubility from Arrhenius if D_0/E_d given |
| `Membrane` | `Membrane` | Metal wall properties; same Arrhenius logic for D, K_S |
| `FluidMaterial` | `FluidMaterial` | Pure property container returned by `Flibe(T)`, `Sodium(T)`, `LiPb(T)` |
| `SolidMaterial` | `SolidMaterial` | Pure property container returned by `Steel(T)` |
| `WireCoil` | `WireCoil` | Wire-coil turbulator with pitch-dependent Sh correlation |
| `CustomTurbulator` | `CustomTurbulator` | User-defined Sh = a·Reᵇ·Scᶜ turbulator |
| `Component` | `Component` | Full PAV/HX permeation component |
| `GLC_Gas` | `GLC_Gas` | Sweep-gas side of a packed column |
| `GLC` | `GLC` | Gas-Liquid Contactor packed column |
| `BreedingBlanket` | `BreedingBlanket` | Tritium-producing breeding blanket |
| `Circuit` | `Circuit` | Series circuit of components |

### Core functions (Component)

| Julia function | Python method | Description |
|---|---|---|
| `use_analytical_efficiency!(c)` | `use_analytical_efficiency()` | Compute eff via Lambert W / closed-form; copy to `c.eff` |
| `analytical_efficiency!(c)` | `analytical_efficiency()` | As above but only fills `c.eff_an` |
| `get_efficiency!(c)` | `get_efficiency()` | Numerical 1-D integration along pipe |
| `get_flux!(c, conc)` | `get_flux(c)` | Permeation flux at local concentration `conc` |
| `outlet_c_comp!(c)` | `outlet_c_comp()` | Fill `c.c_out`, accounting for recirculation |
| `get_regime(c)` | `get_regime()` | Return regime string; optionally print H, W |
| `get_adimensionals!(c)` | `get_adimensionals()` | Compute and store H, W |
| `get_inventory!(c)` | `get_inventory()` | Total tritium inventory (fluid + membrane) |
| `get_solid_inventory!(c)` | `get_solid_inventory()` | Membrane inventory only |
| `get_fluid_inventory!(c)` | `get_fluid_inventory()` | Fluid inventory only |
| `get_pressure_drop!(c)` | `get_pressure_drop()` | Darcy-Weisbach pressure drop |
| `get_pumping_power!(c)` | `get_pumping_power()` | Pumping power [W] |
| `get_global_HX_coeff!(c)` | `get_global_HX_coeff()` | Overall heat transfer coefficient U |
| `inspect(obj)` | `inspect()` | Print all or named fields recursively |
| `update_attribute!(obj, name, val)` | `update_attribute()` | Recursively set a field by name string |

### Material functions

| Julia call | Returns | Fluid type |
|---|---|---|
| `Flibe(T)` | `FluidMaterial` | Molten salt (MS=true) |
| `Sodium(T)` | `FluidMaterial` | Liquid metal (MS=false) |
| `LiPb(T)` | `FluidMaterial` | Liquid metal (MS=false) — mostly TODO placeholders |
| `Steel(T)` | `SolidMaterial` | — |

---

## Differences from the Python version

| Topic | Python | Julia |
|---|---|---|
| Mutation style | Methods mutate `self` in-place | Functions with `!` suffix mutate their first argument |
| Type system | Duck-typed; `isinstance` checks | Abstract type `TriomaClass`; concrete `mutable struct` |
| Pattern matching | `match` / `case` statements | `if / elseif` chains (equivalent semantics) |
| Plotting | `matplotlib` via `plot_component()`, `plot_circuit()` | Not yet ported; plotting methods are omitted |
| `connect_to_component` | Monkey-patched onto classes at import time | Implemented as a standalone `connect_to_component!` function (see `ComponentTools.jl`) |
| `split_HX` / `converge_split_HX` | Returns a `Circuit` | Signature identical; returns a `Circuit` — verify it is present in your build |
| Scipy `integrate.nquad` | Used for 2-D inventory integrals | Replaced by nested `quadgk` calls |
| Scipy `minimize` (Powell) | Used for flux/c_out solvers | Replaced by `Optim.optimize` with `Fminbox(NelderMead())` — worth verifying solver tolerances match for your use case |
| Lambert W | `scipy.special.lambertw` | `LambertW.lambertw` from the LambertW.jl package |

---

## Known limitations and TODOs

- **Plotting** (`plot_component`, `plot_circuit`, `split_HX` convergence plots) is not ported. The Python originals use `matplotlib`; consider `Plots.jl` or `Makie.jl` if you need visualizations.
- **LiPb, Sodium** thermal conductivity and specific heat are placeholder values (`1.0`) inherited from the Python source — they are marked `# TODO` there as well.
- **Steel** Sievert constant `K_S = 1.0` is a placeholder in both the Python and Julia versions.
- **`split_HX` and `converge_split_HX`** are defined in `PAV.jl` but require `Circuit.jl` to be present and loaded.
- The **`BreedingBlanket`** and **`Circuit`** modules follow the same pattern as the other modules but need to be confirmed present in your local build before running circuit-level examples.

---

## Citations

This Julia port implements the physics from the following papers (cited in the original TRIOMA repository):

- Alberghi, C., et al. "Development of new analytical tools for tritium transport modelling." *Fusion Engineering and Design* 177 (2022): 113083.
- Humrickhouse, P. W., and T. F. Fuerst. *Tritium transport phenomena in molten-salt reactors.* INL/EXT-20-59927-Rev000, Idaho National Lab., 2020.
- Fuerst, T. F., C. N. Taylor, and P. W. Humrickhouse. *Molten Salt Tritium Transport Experiment Design.* INL/EXT-21-63108-Rev000, Idaho National Lab., 2021.
- Rader, J. D., M. S. Greenwood, and P. W. Humrickhouse. "Verification of modelica-based models with analytical solutions for tritium diffusion." *Nuclear Technology* 203.1 (2018): 58–65.
- Urgorri, F. R., et al. "Theoretical evaluation of the tritium extraction from liquid metal flows through a free surface and through a permeable membrane." *Nuclear Fusion* (2023).

## License

MIT — same as the upstream TRIOMA Python package.
