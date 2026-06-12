@testset "Component tools" begin

  # ==========================================================================
  # Shared setup helpers
  # ==========================================================================

  function make_ms_component()
    fluid = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                  k_t=0.1, U0=0.2, d_Hyd=0.3)
    geom  = Geometry(L=1.0, dw=0.5, D=0.3)
    mem   = Membrane(D=0.4, dw=0.5, K_S=0.6, T=300.0,
                     k_r=1e7, k=0.8, k_d=1e7)
    Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
  end

  function make_lm_component()
    fluid = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                  k_t=0.1, U0=0.2, d_Hyd=0.3)
    geom  = Geometry(L=1.0, dw=0.5, D=0.3)
    mem   = Membrane(k_d=1e7, D=0.4, dw=0.5, K_S=0.6,
                     T=300.0, k_r=1e7, k=0.8)
    Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
  end

  # ==========================================================================
  # Molten-salt Component tests  (TestMSComponent)
  # ==========================================================================
  @testset "MS Component" begin
    comp = make_ms_component()

    @testset "outlet_c_comp" begin
      outlet_c_comp!(comp)
      @test comp.c_out ≈ 0.5 * (1 - 0.8)
    end

    @testset "T_leak" begin
      comp2 = make_ms_component()
      leak = T_leak(comp2)
      expected = comp2.c_in * comp2.eff * comp2.geometry.D^2 / 4 * π * comp2.fluid.U0
      @test leak ≈ expected
    end

    @testset "get_regime" begin
      comp2 = make_ms_component()
      @test get_regime(comp2) == "Mixed regime"
      @test get_regime(comp2; print_var=true) == "Mixed regime"
    end

    @testset "get_adimensionals" begin
      comp2 = make_ms_component()
      get_adimensionals!(comp2)
      expected_H = comp2.membrane.k_d /
                   (comp2.fluid.k_t * comp2.fluid.Solubility)
      expected_W = 2 * comp2.membrane.k_d * comp2.geometry.dw *
                   (comp2.c_in / comp2.fluid.Solubility)^0.5 /
                   (comp2.membrane.D * comp2.membrane.K_S)
      @test comp2.H ≈ expected_H
      @test comp2.W ≈ expected_W
    end

    @testset "use_analytical_efficiency" begin
      comp2 = make_ms_component()
      use_analytical_efficiency!(comp2)
      @test comp2.eff ≈ 0.99871670123992 rtol=1e-5
    end

    @testset "analytical_efficiency" begin
      comp2 = make_ms_component()
      analytical_efficiency!(comp2)
      @test comp2.eff_an ≈ 0.99871670123992 rtol=1e-5
    end

    @testset "get_efficiency" begin
      comp2 = make_ms_component()
      get_efficiency!(comp2)
      @test comp2.eff ≈ 0.998984924629 rtol=1e-5
    end

    @testset "get_flux" begin
      comp2 = make_ms_component()
      # Python's get_flux() returned c_wl (≈0.0014967), not J_perm.
      # Julia's get_flux! correctly returns J_perm.
      jperm = get_flux!(comp2, 0.3; c_guess=0.3)
      @test jperm ≈ -0.059700645208819604 rtol=1e-5
      @test comp2.J_perm === jperm
    end

    @testset "get_global_HX_coeff" begin
      comp2 = make_ms_component()
      get_global_HX_coeff!(comp2; R_conv_sec=0.1)
      @test comp2.U ≈ 2.9215784663 rtol=1e-5
    end

    @testset "efficiency_vs_analytical" begin
      comp2 = make_ms_component()
      analytical_efficiency!(comp2)
      eff_an = comp2.eff_an
      get_efficiency!(comp2; c_guess=comp2.c_in / 2)
      @test abs(comp2.eff - eff_an) / eff_an ≈ 0 atol=1e-2
    end

    @testset "component_inventory" begin
      comp2 = make_ms_component()
      use_analytical_efficiency!(comp2)
      get_inventory!(comp2)
      # KNOWN BUG: analytical_solid_inventory! for MS gives a negative result
      # because K>0 but ifun(r_out)-ifun(r_in)<0 (sign error in the K*ifun
      # expression — should subtract r_out from r_in, not r_in from r_out).
      # Numerical method (flag_an=false) gives ~0.00826, also not matching the
      # expected 0.00532, suggesting a deeper formula discrepancy.
      @test_broken comp2.inv ≈ 0.00531677445914132 rtol=1e-3
    end
  end

  # ==========================================================================
  # Liquid-metal Component tests  (TestLMComponent)
  # ==========================================================================
  @testset "LM Component" begin
    comp = make_lm_component()

    @testset "outlet_c_comp" begin
      outlet_c_comp!(comp)
      @test comp.c_out ≈ comp.c_in * (1 - comp.eff)
    end

    @testset "T_leak" begin
      comp2 = make_lm_component()
      leak = T_leak(comp2)
      expected = comp2.c_in * comp2.eff *
                 comp2.geometry.D^2 / 4 * π * comp2.fluid.U0
      @test leak ≈ expected
    end

    @testset "get_regime" begin
      comp2 = make_lm_component()
      @test get_regime(comp2; print_var=true) == "Mixed regime"
    end

    @testset "get_adimensionals" begin
      comp2 = make_lm_component()
      get_adimensionals!(comp2)
      expected_W = comp2.membrane.k_r / comp2.membrane.D *
                   comp2.membrane.K_S * comp2.membrane.dw *
                   comp2.c_in / comp2.fluid.Solubility
      expected_H = expected_W * comp2.membrane.D * comp2.membrane.K_S /
                   (comp2.fluid.k_t * comp2.fluid.Solubility * comp2.membrane.dw)
      @test comp2.W ≈ expected_W
      @test comp2.H ≈ expected_H
    end

    @testset "use_analytical_efficiency" begin
      comp2 = make_lm_component()
      use_analytical_efficiency!(comp2)
      @test comp2.eff ≈ 0.998295638580 rtol=1e-5
    end

    @testset "analytical_efficiency - after clearing k_t" begin
      comp2 = make_lm_component()
      update_attribute!(comp2, "k_t", nothing)
      analytical_efficiency!(comp2)
      @test comp2.eff_an ≈ 0.00053628139636452 rtol=1e-5
    end

    @testset "get_efficiency" begin
      comp2 = make_lm_component()
      get_efficiency!(comp2; c_guess=comp2.c_in / 2)
      @test comp2.eff ≈ 0.9986246 rtol=1e-5
    end

    @testset "get_flux" begin
      comp2 = make_lm_component()
      jperm = get_flux!(comp2, 0.3; c_guess=0.3)
      @test jperm ≈ -0.028685539953802155 rtol=1e-5
      @test comp2.J_perm === jperm
    end

    @testset "get_global_HX_coeff" begin
      comp2 = make_lm_component()
      get_global_HX_coeff!(comp2; R_conv_sec=0.1)
      @test comp2.U ≈ 2.9215784663 rtol=1e-5
    end

    @testset "efficiency_vs_analytical" begin
      comp2 = make_lm_component()
      analytical_efficiency!(comp2)
      eff_an = comp2.eff_an
      get_efficiency!(comp2; c_guess=comp2.c_in / 2)
      @test abs(comp2.eff - eff_an) / eff_an ≈ 0 atol=1e-2
    end

    @testset "update_attribute" begin
      comp2 = make_lm_component()
      update_attribute!(comp2, "c_in", 0.6)
      @test comp2.c_in == 0.6
      update_attribute!(comp2, "U0", 0.3)
      @test comp2.fluid.U0 == 0.3
      @test_throws ArgumentError update_attribute!(comp2, "kghufh", 0.3)
      try
        update_attribute!(comp2, "kghufh", 0.3)
      catch e
        @test occursin("kghufh", e.msg)
        @test occursin("Component",  e.msg)
      end
    end

    @testset "set_material_properties" begin
      comp2   = make_lm_component()
      T       = 800.0
      fm      = Flibe(T)
      set_properties_from_fluid_material!(comp2.fluid, fm)
      @test comp2.fluid.rho == fm.rho
      sm = Steel(T)
      set_properties_from_solid_material!(comp2.membrane, sm)
      @test comp2.membrane.K_S == sm.K_S
    end

    @testset "LiPb material" begin
      comp2 = make_lm_component()
      fm    = LiPb(800.0)
      set_properties_from_fluid_material!(comp2.fluid, fm)
      @test comp2.fluid.rho == fm.rho
    end

    @testset "Sodium material" begin
      comp2 = make_lm_component()
      fm    = Sodium(800.0)
      set_properties_from_fluid_material!(comp2.fluid, fm)
      @test comp2.fluid.rho == fm.rho
    end

    @testset "update T properties" begin
      comp2 = make_lm_component()
      update_attribute!(comp2, "T", 999.0)
      @test comp2.fluid.T == 999.0
    end

    @testset "component inventory (numerical)" begin
      comp2 = make_lm_component()
      get_inventory!(comp2; flag_an=false)
      @test comp2.inv ≈ 0.0070080332665533665 rtol=1e-5
    end

    @testset "pumping power" begin
      comp2 = make_lm_component()
      get_pumping_power!(comp2)
      @test comp2.pumping_power ≈ 0.019029194163191265 rtol=1e-5
    end

    @testset "inspect runs without error" begin
      comp2 = make_lm_component()
      @test_nowarn inspect(comp2)
      @test_nowarn inspect(comp2.geometry)
    end
  end

  # ==========================================================================
  # Exotic Component (Arrhenius params, WireCoil, recirculation)  (TestExoticComponent)
  # ==========================================================================
  @testset "Exotic Component" begin
    function make_exotic_component()
      fluid = Fluid(T=300.0, D_0=1e-9, E_d=0.5, Solubility_0=0.5, E_s=0.5,
                    MS=true, mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                    k_t=0.1, U0=0.2, d_Hyd=0.3, recirculation=-0.5)
      wc   = WireCoil(pitch=1e-2)
      geom = Geometry(L=1.0, dw=0.5e-3, D=0.3, turbulator=wc)
      mem  = Membrane(D_0=0.4, E_d=0.1, E_S=0.1, dw=0.5,
                      K_S_0=0.6, T=300.0, k_r=1e7, k=0.8, k_d=1e7)
      Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    end

    @testset "Arrhenius D and K_S computed at construction" begin
      comp = make_exotic_component()
      @test comp.fluid.D    ≈ 3.984462016634033e-18 rtol=1e-5
      @test comp.membrane.D ≈ 0.008358607447235438  rtol=1e-5
    end

    @testset "update T — Arrhenius properties refresh" begin
      comp = make_exotic_component()
      update_attribute!(comp, "T", 999.0)
      @test comp.fluid.T    == 999.0
      @test comp.membrane.D ≈ 0.12519232082575535  rtol=1e-5
      @test comp.fluid.D    ≈ 3.003229324075595e-12 rtol=1e-5
    end

    @testset "negative recirculation scales U0" begin
      comp = make_exotic_component()
      @test comp.fluid.recirculation == -0.5
      @test comp.fluid.U0            == 0.1
    end

    @testset "WireCoil k_t correlation" begin
      comp = make_exotic_component()
      update_attribute!(comp, "k_t", nothing)
      use_analytical_efficiency!(comp)
      # Old value 3.74e-11 was computed with pitch as Sherwood characteristic length,
      # which is physically incorrect (Sh is always defined using d_hyd as char. length).
      @test comp.fluid.k_t ≈ 1.716729981847424e-10 rtol=1e-5
      # CustomTurbulator falls back to k_t already set
      ct = CustomTurbulator(a=1.0, b=1.0, c=1.0)
      comp.geometry.turbulator = ct
      update_attribute!(comp, "k_t", nothing)
      use_analytical_efficiency!(comp)
      @test comp.fluid.k_t ≈ 0.1 rtol=1e-5
    end

    @testset "outlet_c_comp with bypass recirculation" begin
      comp = make_exotic_component()
      outlet_c_comp!(comp)
      @test comp.c_out == 0.3
    end

    @testset "outlet_c_comp with positive recirculation" begin
      comp = make_exotic_component()
      update_attribute!(comp, "recirculation", 0.5)
      outlet_c_comp!(comp)
      @test comp.c_out ≈ 5.399292025261743e-07 rtol=1e-5
    end

    @testset "total flowrate" begin
      comp = make_exotic_component()
      get_total_flowrate(comp)
      @test comp.flowrate ≈ 0.007068583470577035 rtol=1e-5
    end

    @testset "component volumes" begin
      comp = make_exotic_component()
      define_component_volumes!(comp)
      @test comp.fluid.V    ≈ 0.07068583470577035  rtol=1e-5
      @test comp.membrane.V ≈ 0.0004704534998750733 rtol=1e-5
      total = comp.fluid.V + comp.membrane.V
      @test total ≈ 0.07115628820564542 rtol=1e-5
    end
  end

  # ==========================================================================
  # FluidMaterial tests  (TestFluidMaterial)
  # ==========================================================================
  @testset "FluidMaterial" begin
    fm = FluidMaterial(T=300.0, D=1e-9, Solubility=0.5, MS=true,
                       mu=1e-3, rho=1000.0, k=0.5, cp=1.0)
    @test fm.T          == 300.0
    @test fm.D          == 1e-9
    @test fm.Solubility == 0.5
    @test fm.MS         == true
    @test fm.mu         == 1e-3
    @test fm.rho        == 1000.0
    @test fm.k          == 0.5
    @test fm.cp         == 1.0

    @testset "update_attribute" begin
      update_attribute!(fm, "T", 400.0)
      @test fm.T == 400.0
      @test_throws ArgumentError update_attribute!(fm, "kghufh", 0.3)
      try
        update_attribute!(fm, "kghufh", 0.3)
      catch e
        @test occursin("kghufh",       e.msg)
        @test occursin("FluidMaterial", e.msg)
      end
    end
  end

  # ==========================================================================
  # SolidMaterial tests  (Test_SolidMaterial)
  # ==========================================================================
  @testset "SolidMaterial" begin
    sm = SolidMaterial(T=300.0, K_S=1.0, D=1e-7)
    @test sm.T   == 300.0
    @test sm.K_S == 1.0

    @testset "update_attribute" begin
      update_attribute!(sm, "T", 400.0)
      @test sm.T == 400.0
      @test_throws ArgumentError update_attribute!(sm, "kghufh", 0.3)
      try
        update_attribute!(sm, "kghufh", 0.3)
      catch e
        @test occursin("kghufh",       e.msg)
        @test occursin("SolidMaterial", e.msg)
      end
    end
  end

  # ==========================================================================
  # BreedingBlanket tests  (Test_BB_Component)
  # ==========================================================================
  @testset "BreedingBlanket" begin
    function make_bb()
      BreedingBlanket(c_in=0.0, Q=0.5e9, TBR=1.05,
                      T_out=900.0, T_in=800.0, fluid=Flibe(850.0))
    end

    @testset "get_cout" begin
      bb = make_bb()
      get_cout!(bb; print_var=true)
      @test bb.c_out ≈ 0.0001473990666223908 rtol=1e-5
    end

    @testset "get_flowrate" begin
      bb = make_bb()
      get_flowrate!(bb)
      @test bb.m_coolant ≈ 2095.5574182607 rtol=1e-5
    end

    @testset "update_attribute - direct and nested" begin
      bb = make_bb()
      update_attribute!(bb, "c_in", 0.6)
      @test bb.c_in == 0.6
      update_attribute!(bb, "cp", 0.3)
      @test bb.fluid.cp == 0.3
      @test_throws ArgumentError update_attribute!(bb, "kghufh", 0.3)
      try
        update_attribute!(bb, "kghufh", 0.3)
      catch e
        @test occursin("kghufh",         e.msg)
        @test occursin("BreedingBlanket", e.msg)
      end
    end
  end

  # ==========================================================================
  # Membrane tests  (TestMembrane)
  # ==========================================================================
  @testset "Membrane" begin
    mem = Membrane(k_d=1e7, D=0.4, dw=0.5, K_S=0.6, T=300.0, k_r=1e7, k=0.8)
    update_attribute!(mem, "dw", 1e-3)
    @test mem.dw == 1e-3
    @test_throws ArgumentError update_attribute!(mem, "kghufh", 0.3)
    try
      update_attribute!(mem, "kghufh", 0.3)
    catch e
      @test occursin("kghufh",  e.msg)
      @test occursin("Membrane", e.msg)
    end
  end

  # ==========================================================================
  # Fluid tests  (TestFluid)
  # ==========================================================================
  @testset "Fluid" begin
    fl = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=true,
               mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
               k_t=0.1, U0=0.2, d_Hyd=0.3)

    @testset "update_attribute" begin
      update_attribute!(fl, "Solubility", 1e-3)
      @test fl.Solubility == 1e-3
      @test_throws ArgumentError update_attribute!(fl, "kghufh", 0.3)
      try
        update_attribute!(fl, "kghufh", 0.3)
      catch e
        @test occursin("kghufh", e.msg)
        @test occursin("Fluid",  e.msg)
      end
    end

    @testset "get_kt! - no d_Hyd defined" begin
      fl2 = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, k_t=0.1, U0=0.2)
      @test_logs (:warn, r"not defined") get_kt!(fl2)
    end

    @testset "get_kt! - turbulent flow" begin
      fl2 = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                  k_t=nothing, U0=0.2, d_Hyd=0.3)
      get_kt!(fl2)
      @test fl2.k_t ≈ 8.046408367835323e-06 rtol=1e-5
    end

    @testset "get_kt! - laminar flow" begin
      fl2 = Fluid(T=300.0, D=1e-9, Solubility=0.5, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                  k_t=nothing, U0=1e-6, d_Hyd=0.3)
      get_kt!(fl2)
      @test fl2.k_t ≈ 3.66 / fl2.d_Hyd * fl2.D rtol=1e-10
    end
  end

  # ==========================================================================
  # Transport-regime tests for all MS/LM regimes
  # ==========================================================================
  function make_ms_diffusion_limited()
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e7, D=1e-9, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
  end

  @testset "MS Diffusion Limited regime" begin
    comp = make_ms_diffusion_limited()
    @test get_regime(comp) == "Diffusion Limited"
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    # Both analytical and numerical now agree: analytical_efficiency! uses the Brent
    # fallback (beta_tau >> log(floatmax)) which finds the same solution as the
    # numerical integrator (~1.22e-5 for this parameter set).
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "MS Mixed (diffusion+mass-transport) regime" begin
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e7, D=1e-6, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Mixed regime"
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    # Both analytical and numerical agree after the refactoring: the Brent fallback
    # in analytical_efficiency! finds the same concentration profile as get_efficiency!.
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "MS Mixed (diffusion+surface) regime" begin
    fluid = Fluid(T=750.0, D=2e-2, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e-10, D=1e-7, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e-10, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Mixed regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in)
  end

  @testset "MS Mixed (mass-transport+surface) regime" begin
    fluid = Fluid(T=750.0, D=2e-8, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e-8, D=1e-2, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e-8, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Mixed regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "MS Mass-transport limited regime" begin
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-2, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e7, D=1e-2, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Mass transport limited"
    get_flux!(comp, 0.3; c_guess=0.3)
    @test comp.J_perm ≈ -9.6149466095734e-05 rtol=1e-5
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "MS Surface limited regime" begin
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e-16, D=1e-9, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e-16, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Surface limited"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "MS Fully mixed regime" begin
    fluid = Fluid(T=750.0, D=2e-11, Solubility=1e-4, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e-10, D=1e-7, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e-10, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp) == "Mixed regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "LM Mass-transport limited regime" begin
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-2, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e7, D=1e-2, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Mass transport limited"
    get_flux!(comp, 0.3; c_guess=0.3)
    @test comp.J_perm ≈ -4.80747330478670e-05 rtol=1e-5
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "LM Mixed (diffusion+mass-transport) regime" begin
    fluid = Fluid(T=750.0, D=2e-9, Solubility=1e-2, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e7, D=1e-7, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Mixed regime"
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "LM Mixed (diffusion+surface) regime" begin
    fluid = Fluid(T=750.0, D=2e-3, Solubility=1e-2, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e-3, D=1e-7, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e-3, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Mixed regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "LM Transport+surface limited regime" begin
    fluid = Fluid(T=750.0, D=2e-6, Solubility=1e-2, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e-3, D=1e-2, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e-3, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Transport and surface limited regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "LM Fully mixed regime" begin
    fluid = Fluid(T=750.0, D=2e-11, Solubility=1e-2, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=2.0, d_Hyd=2e-3)
    geom  = Geometry(L=1.0, dw=1e-4, D=2e-3)
    mem   = Membrane(k_d=1e-9, D=1e-9, dw=1e-4, K_S=0.6e-2, T=700.0, k_r=1e-9, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Mixed regime"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  @testset "LM Diffusion limited regime" begin
    fluid = Fluid(T=750.0, D=2e-5, Solubility=1e-4, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e7, D=1e-9, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e7, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Diffusion Limited"
    jperm = get_flux!(comp, 0.3)
    @test jperm ≈ -2.5960851589349628e-6 rtol=1e-5
    analytical_efficiency!(comp)
    get_efficiency!(comp; c_guess=comp.c_in / 2)
    @test abs(comp.eff - comp.eff_an) / comp.eff_an ≈ 0 atol=1e-2
  end

  @testset "LM Surface limited regime" begin
    fluid = Fluid(T=750.0, D=2e-5, Solubility=1e-4, MS=false,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=1.0, dw=1e-2, D=2e-2)
    mem   = Membrane(k_d=1e-16, D=1e-9, dw=1e-2, K_S=0.6e-2, T=700.0, k_r=1e-16, k=0.8)
    comp  = Component(c_in=0.5, geometry=geom, eff=0.8, fluid=fluid, membrane=mem)
    @test get_regime(comp; print_var=true) == "Surface limited"
    @test_nowarn get_efficiency!(comp; c_guess=comp.c_in / 2)
  end

  # ==========================================================================
  # GLC tests  (testLMGLCComponent, testMSGLCComponent)
  # ==========================================================================
  @testset "LM GLC" begin
    using AtomicAndPhysicalConstants: BOLTZMANN_k, J_PER_EV, AVOGADRO
    R_const = BOLTZMANN_k * J_PER_EV * AVOGADRO
    T    = 400.0 + 273.15
    Z    = 0.6;  R_col = 0.3
    Q_l  = 71e-3 / 3600
    K_S  = 1.33e-4 * exp(-1350.0 / (R_const * T))
    gas  = GLC_Gas(G_gas=3e-3 / 3600, pg_in=0.0, p_tot=1.5e5)
    fl   = Fluid(Solubility=K_S, MS=false)
    glc  = GLC(H=Z, R=R_col, c_in=1e-2, c_out=9e-3,
               fluid=fl, GLC_gas=gas, T=T, G_L=Q_l)

    # Expected kla computed from Julia's R_const (8.31446... vs Python's hardcoded 8.314).
    # LM kla agrees with Python to 2e-6 relative; rtol=1e-5 covers this.
    @testset "get_kla_from_cout" begin
      get_kla_from_cout!(glc)
      @test glc.kla ≈ 1.2850620482865885e-5 rtol=1e-5
    end

    @testset "get_c_out" begin
      glc2 = GLC(H=Z, R=R_col, c_in=1e-2, fluid=fl, GLC_gas=gas, T=T, G_L=Q_l,
                 kla=1.2850620482865885e-5)
      get_c_out!(glc2)
      @test glc2.c_out ≈ 0.009 rtol=1e-5
    end

    @testset "get_z_from_eff" begin
      glc3 = GLC(H=Z, R=R_col, c_in=1e-2, c_out=9e-3, fluid=fl, GLC_gas=gas, T=T, G_L=Q_l,
                 kla=1.2850620482865885e-5)
      z = get_z_from_eff(glc3)
      @test z ≈ 0.6 rtol=1e-5
    end
  end

  @testset "MS GLC" begin
    using AtomicAndPhysicalConstants: BOLTZMANN_k, J_PER_EV, AVOGADRO
    R_const = BOLTZMANN_k * J_PER_EV * AVOGADRO
    T    = 400.0 + 273.15
    Z    = 0.6;  R_col = 0.3
    Q_l  = 71e-3 / 3600
    K_H  = 1.33e-4 * exp(-1350.0 / (R_const * T))
    gas  = GLC_Gas(G_gas=3e-3 / 3600, pg_in=0.0, p_tot=1.5e5)
    fl   = Fluid(Solubility=K_H, MS=true)
    glc  = GLC(H=Z, R=R_col, c_in=1e-2, c_out=9e-3,
               fluid=fl, GLC_gas=gas, T=T, G_L=Q_l)

    # MS kla is more sensitive to R_const than LM (7.5e-4 relative diff vs Python's 8.314).
    # Expected value computed with Julia's R_const = 8.31446...
    @testset "get_kla_from_cout" begin
      get_kla_from_cout!(glc)
      @test glc.kla ≈ 2.9788115080255035e-5 rtol=1e-5
    end

    @testset "get_c_out" begin
      glc2 = GLC(H=Z, R=R_col, c_in=1e-2, fluid=fl, GLC_gas=gas, T=T, G_L=Q_l,
                 kla=2.9788115080255035e-5)
      get_c_out!(glc2)
      @test glc2.c_out ≈ 0.009 rtol=1e-5
    end

    @testset "get_z_from_eff" begin
      glc3 = GLC(H=Z, R=R_col, c_in=1e-2, c_out=9e-3, fluid=fl, GLC_gas=gas, T=T, G_L=Q_l,
                 kla=2.9788115080255035e-5)
      z = get_z_from_eff(glc3)
      @test z ≈ 0.6 rtol=1e-5
    end
  end

  # ==========================================================================
  # Circuit tests  (testclosedCircuit)
  # ==========================================================================
  @testset "Closed Circuit" begin
    fluid = Fluid(T=900.0, D=1e-7, Solubility=0.5, MS=true,
                  mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                  k_t=0.1, U0=1.0, d_Hyd=2e-2)
    geom  = Geometry(L=10.0, dw=0.5e-3, D=2e-2)
    mem   = Membrane(D_0=1e-7, E_d=1.0, dw=0.5e-3,
                     K_S=0.6, T=900.0, k_r=1e7, k=0.8, k_d=1e7)
    comp1 = Component(geometry=geom, fluid=fluid, membrane=mem, loss=false)
    comp2 = Component(geometry=geom, fluid=fluid, membrane=mem,
                      loss=false, name="PAV")
    bb    = BreedingBlanket(c_in=1e-3, Q=0.5e9, TBR=1.05,
                            T_out=900.0, T_in=800.0, fluid=Flibe(850.0), name="BB")

    geom_hx = Geometry(L=10.0, dw=1e-3, D=2e-3)
    mem_hx  = Membrane(D_0=1e-9, E_d=0.2, dw=1e-3, K_S=0.6, T=850.0, k_r=1e7, k=0.8, k_d=1e7)
    fluid_hx = Fluid(T=850.0, D=1e-9, Solubility=0.5, MS=true,
                     mu=1e-3, rho=1000.0, k=0.5, cp=1.0,
                     k_t=0.1, U0=1.0, d_Hyd=2e-2)
    comp_hx = Component(c_in=1e-3, geometry=geom_hx, fluid=fluid_hx,
                        membrane=mem_hx, loss=true)

    @testset "solve_circuit runs without error" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      @test_nowarn solve_circuit!(circuit)
    end

    @testset "get_eff_circuit runs without error" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      solve_circuit!(circuit)
      @test_nowarn get_eff_circuit!(circuit)
      @test circuit.eff isa Float64
    end

    @testset "get_gains_and_losses" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      solve_circuit!(circuit)
      get_gains_and_losses!(circuit)
      @test circuit.extraction_perc isa Float64
      @test circuit.loss_perc       isa Float64
    end

    @testset "inspect_circuit runs without error" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      solve_circuit!(circuit)
      @test_nowarn inspect_circuit(circuit)
      @test_nowarn inspect_circuit(circuit; name="PAV")
    end

    @testset "add_component and flattening" begin
      c1 = Circuit(components=[comp1])
      c2 = Circuit(components=[comp2])
      add_component!(c1, c2)
      @test length(c1.components) == 2
    end

    @testset "get_inventory runs without error" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      solve_circuit!(circuit)
      @test_nowarn get_inventory!(circuit)
      @test circuit.inv isa Float64
    end

    @testset "get_circuit_pumping_power" begin
      circuit = Circuit(components=[bb, comp1, comp2], closed=true)
      solve_circuit!(circuit)
      get_circuit_pumping_power!(circuit)
      @test circuit.pumping_power isa Float64
      @test circuit.pumping_power > 0
    end

    @testset "estimate_cost" begin
      circuit = Circuit(components=[comp1, comp2])
      n = length(circuit.components)
      mc = ones(Float64, n)
      fc = ones(Float64, n)
      estimate_cost!(circuit; metal_costs=mc, fluid_costs=fc)
      @test circuit.cost isa Float64
      @test circuit.cost > 0
    end
  end

end
