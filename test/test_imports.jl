@testset "Imports" begin

  # All of these should be accessible after `using TRIOMA`

  @testset "Core types importable" begin
    @test TriomaClass       isa DataType
    @test Geometry          isa DataType
    @test Fluid             isa DataType
    @test Membrane          isa DataType
    @test FluidMaterial     isa DataType
    @test SolidMaterial     isa DataType
    @test WireCoil          isa DataType
    @test CustomTurbulator  isa DataType
    @test Component         isa DataType
    @test GLC_Gas           isa DataType
    @test GLC               isa DataType
    @test BreedingBlanket   isa DataType
    @test Circuit           isa DataType
  end

  @testset "Material functions importable" begin
    @test Flibe   isa Function
    @test Sodium  isa Function
    @test LiPb    isa Function
    @test Steel   isa Function
  end

  @testset "Component functions importable" begin
    @test use_analytical_efficiency!   isa Function
    @test analytical_efficiency!       isa Function
    @test get_efficiency!              isa Function
    @test outlet_c_comp!               isa Function
    @test get_regime                   isa Function
    @test get_inventory!               isa Function
    @test get_pressure_drop!           isa Function
    @test get_pumping_power!           isa Function
    @test get_global_HX_coeff!        isa Function
  end

  @testset "Circuit functions importable" begin
    @test solve_circuit!               isa Function
    @test get_eff_circuit!             isa Function
    @test get_gains_and_losses!        isa Function
    @test get_circuit_pumping_power!   isa Function
    @test add_component!               isa Function
    @test connect_to_component!        isa Function
    @test inspect_circuit              isa Function
  end

  @testset "Correlation functions importable" begin
    @test Nu_DittusBoelter   isa Function
    @test Sherwood           isa Function
    @test Re                 isa Function
    @test Pr                 isa Function
    @test Schmidt            isa Function
    @test get_deltaTML       isa Function
  end

  @testset "BreedingBlanket functions importable" begin
    @test get_flowrate!  isa Function
    @test get_cout!      isa Function
  end

  @testset "GLC functions importable" begin
    @test get_c_out!          isa Function
    @test get_kla_from_cout!  isa Function
    @test get_z_from_eff      isa Function
  end

end
