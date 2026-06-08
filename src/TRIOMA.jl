"""
    TRIOMA

Julia port of the TRIOMA (TRItium Object-oriented and Modular Analysis) Python package.

Original Python package by Gabriele Ferrero & Samuele Meschini (Politecnico di Torino).
MIT License.

## Package structure

| Julia module           | Python source file            |
|------------------------|-------------------------------|
| `Correlations`         | `correlations.py`             |
| `PipeSubclasses`       | `PipeSubclasses.py`|
| `MoltenSalts`          | `molten_salts.py`             |
| `LiquidMetals`         | `liquid_metals.py`            |
| `Materials`            | `materials.py`                |
| `ExtractorFunctions`   | `extractor.py`     |
| `GasLiquidContactor`   | `GasLiquidContactor.py` |
| `PAVComponent`         | `PAV.py`           |
| `BreedingBlanketModule`| `BreedingBlanket.py`          |
| `CircuitModule`        | `Circuit.py`                  |

## Key types

- `Component`       – PAV / heat-exchanger pipe
- `BreedingBlanket` – tritium-generating blanket component
- `GLC`             – Gas-Liquid Contactor (packed column)
- `Circuit`         – series chain of components
- `Fluid`, `Membrane`, `Geometry` – sub-objects attached to `Component`
- `FluidMaterial`, `SolidMaterial` – property containers
- `WireCoil`, `CustomTurbulator`   – turbulator types

## Notes on translation

1. Python uses mutable class instances; Julia uses `mutable struct`.
2. Python `None` → Julia `nothing`; all nullable fields are `Union{T,Nothing}`.
3. Python `match ... case` → Julia `if/elseif` chains.
4. `scipy.optimize.minimize(method="Powell")` → `Optim.jl` with `Fminbox(Powell())`.
5. `scipy.integrate.quad / fixed_quad / nquad` → `QuadGK.jl`.
6. `scipy.special.lambertw` → `SpecialFunctions.lambertw`.
7. Plotting (matplotlib) is **not** translated; visualisation is left to the user.
8. The `split_HX` / `converge_split_HX` methods are complex and depend on
   `Circuit`; they are omitted here.  Worth implementing separately once the
   rest of the package is validated.

## Usage example

```julia
using TRIOMA
using TRIOMA.PipeSubclasses, TRIOMA.PAVComponent

# Build a simple PAV component
geom = Geometry(L=1.0, D=0.01, thick=5e-4)
fluid = Fluid(T=700.0, D=1e-9, Solubility=1e-3, MS=true, mu=1e-3,
              rho=2000.0, U0=0.5, d_Hyd=0.01, k=1.1, cp=2386.0)
mem   = Membrane(T=700.0, D=1e-12, thick=5e-4, K_S=1e5, k_d=1e-4)
comp  = Component(geometry=geom, fluid=fluid, membrane=mem, c_in=1e-3)

# Evaluate analytical efficiency
analytical_efficiency!(comp)
println("Efficiency: ", comp.eff_an)
```
"""
module TRIOMA

# ── Re-export the abstract base and its helpers under a stable module name ─────
module TriomaModule
    include("TriomaClass.jl")
    export TriomaClass, inspect, update_attribute!
end

include("Correlations.jl")

# Make sub-modules visible to each other via parent module
using .TriomaModule

include("PipeSubclasses.jl")
# PipeSubclasses references TriomaModule and Correlations via parent-module lookup.
# We give it explicit access by re-exporting.

include("RegimeHelpers.jl")
include("Materials.jl")
include("ExtractorFunctions.jl")
include("GasLiquidContactor.jl")
include("PAVComponent.jl")
include("BreedingBlanket.jl")
include("Circuit.jl")

# ── Public re-exports ──────────────────────────────────────────────────────────
using .TriomaModule: TriomaClass, inspect, update_attribute!
using .Correlations
using .PipeSubclasses
using .MoltenSalts
using .LiquidMetals
using .Materials
using .ExtractorFunctions
using .GasLiquidContactor
using .PAVComponent
using .BreedingBlanketModule
using .CircuitModule

export
    # Base
    TriomaClass, inspect, update_attribute!,
    # Sub-types
    Geometry, Fluid, Membrane, FluidMaterial, SolidMaterial,
    Turbulator, WireCoil, CustomTurbulator,
    Component, BreedingBlanket, GLC, GLC_Gas, Circuit,
    # Material factories
    Flibe, Sodium, LiPb, Steel,
    # Key methods (bang = mutates)
    get_kt!, update_T_prop!,
    analytical_efficiency!, use_analytical_efficiency!,
    get_efficiency!, get_flux!, outlet_c_comp!,
    get_regime, get_adimensionals!,
    get_pressure_drop!, get_pumping_power!,
    get_pipe_flowrate, get_total_flowrate,
    get_solid_inventory!, get_fluid_inventory!, get_inventory!,
    get_global_HX_coeff!, T_leak,
    get_flowrate!, get_cout!,
    get_c_out!, get_kla_from_cout!, get_z_from_eff,
    get_eff_circuit!, solve_circuit!, get_circuit_inventory!,
    get_circuit_pumping_power!, get_gains_and_losses!,
    add_component!, inspect_circuit,
    # Correlations
    Nu_SiederTate, Nu_Gnielinsky, Nu_DittusBoelter, f_Pethukov, get_h_from_Nu,
    f_Haaland, Schmidt, Sherwood, get_k_from_Sh, Re, Pr,
    Sherwood_bubbles, get_length_HX, get_deltaTML,
    # Regime helpers
    W_ms, H_ms, get_regime_ms, W_lm, partition_param_lm, get_regime_lm,
    # Extractor functions
    NTU_lm, NTU_ms, extractor_lm, extractor_ms,
    length_extractor_lm, length_extractor_ms,
    get_c_out_GLC_lm, get_c_out_GLC_ms,
    calculate_gas_velocity, pack_corr

end # module TRIOMA
