"""
TRIOMA.jl

TRItium Object-oriented and Modular Analysis — Julia port.

A tool for designing the Outer Fuel Cycle (OFC) of fusion reactors.
Provides object-oriented models for tritium transport analysis including
Permeation Against Vacuum (PAV) extractors, Gas-Liquid Contactors,
Breeding Blankets, Heat Exchangers, and full circuit analysis.

Original Python package by Gabriele Ferrero and Samuele Meschini
(Politecnico di Torino). MIT licence.

Usage:
    using TRIOMA

All public types and functions are exported from this module. See the
individual submodule files for full documentation.
"""
module TRIOMA

# ---------------------------------------------------------------------------
# Load submodules in dependency order
# ---------------------------------------------------------------------------

include("TriomaTypes.jl")
include("Correlations.jl")
include("PipeSubclasses.jl")
include("Materials.jl")
include("ExtractorFunctions.jl")
include("RegimeHelpers.jl")
include("GasLiquidContactor.jl")
include("PAV.jl")
include("BreedingBlanket.jl")
include("Circuit.jl")

# ---------------------------------------------------------------------------
# Bring submodules into scope
# ---------------------------------------------------------------------------

using .TriomaTypes
using .Correlations
using .PipeSubclasses
using .Materials
using .FusionCoolant
using .ExtractorFunctions
using .GasLiquidContactor
using .PAVPipe
using .BreedingBlanketModule
using .CircuitModule

# ---------------------------------------------------------------------------
# Re-export the full public API
# ---------------------------------------------------------------------------

# --- Base infrastructure ---
export inspect, update_attribute!

# --- Correlations ---
export Nu_SiederTate, Nu_Gnielinsky, Nu_DittusBoelter, f_Pethukov,
       get_h_from_Nu, f_Haaland, Schmidt, Sherwood, Sherwood_HT_analogy,
       get_k_from_Sh, Re, Pr, Sherwood_bubbles, get_length_HX, get_deltaTML

# --- Pipe subclasses (types) ---
export Geometry, Fluid, Membrane, FluidMaterial, SolidMaterial,
       Turbulator, WireCoil, CustomTurbulator,
       get_fluid_volume, get_solid_volume, get_total_volume

# --- Pipe subclasses (functions) ---
export set_properties_from_fluid_material!, set_properties_from_solid_material!,
       update_T_prop!, get_kt!, k_t_correlation, h_t_correlation

# --- Materials ---
export Flibe, Sodium, LiPb, Steel

# --- Coolant regime dimensionless numbers (molten salt + liquid metal) ---
export W_ms, H_ms, get_regime_ms,
       W_lm, partition_param_lm, get_regime_lm

# --- Packed-column / GLC functions ---
export calculate_gas_velocity,
       NTU_lm, extractor_lm, length_extractor_lm,
       NTU_ms, extractor_ms, length_extractor_ms,
       get_c_out_GLC_lm, get_c_out_GLC_ms,
       pack_corr, corr_packed

# --- GLC types and methods ---
export GLC_Gas, GLC, get_c_out!, get_kla_from_cout!, get_z_from_eff

# --- Component (PAV / HX pipe) ---
export Component,
       get_regime,
       get_pipe_flowrate, get_total_flowrate,
       define_component_volumes!,
       get_adimensionals!,
       analytical_efficiency!, use_analytical_efficiency!,
       get_efficiency!, get_flux!,
       outlet_c_comp!,
       get_global_HX_coeff!,
       get_pressure_drop!, get_pumping_power!,
       get_solid_inventory!, get_fluid_inventory!, get_inventory!,
       analytical_solid_inventory!, analytical_fluid_inventory!,
       estimate_cost!,
       T_leak,
       friction_factor

# --- BreedingBlanket ---
export BreedingBlanket, get_flowrate!, get_cout!, connect_to_component!

# --- Circuit ---
export Circuit,
       add_component!,
       solve_circuit!,
       get_eff_circuit!,
       get_gains_and_losses!,
       get_circuit_pumping_power!,
       inspect_circuit

end # module TRIOMA
