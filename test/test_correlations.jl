@testset "Correlations" begin

  @testset "get_deltaTML" begin
    T_in_hot   = 100.0; T_out_hot  = 80.0
    T_in_cold  = 40.0;  T_out_cold = 20.0
    expected = ((T_in_hot - T_out_cold) - (T_out_hot - T_in_cold)) /
               log((T_in_hot - T_out_cold) / (T_out_hot - T_in_cold))
    @test get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold) == expected
  end

  @testset "Nu_SiederTate" begin
    Re_v = 10000.0; Pr_v = 0.7; mu_w = 0.001; mu_c = 0.002
    expected = 0.027 * Re_v^(4/5) * Pr_v^0.3 * (mu_c / mu_w)^0.14
    @test Nu_SiederTate(Re_v, Pr_v, mu_w, mu_c) == expected
  end

  @testset "Nu_Gnielinsky" begin
    Re_v = 10000.0; Pr_v = 0.7; f = 0.02
    expected = (f/8) * (Re_v - 1000) * Pr_v /
               (1 + 12.7 * (f/8)^0.5 * (Pr_v^(2/3) - 1))
    @test Nu_Gnielinsky(Re_v, Pr_v, f) == expected
  end

  @testset "Nu_DittusBoelter" begin
    Re_v = 10000.0; Pr_v = 0.7
    expected = 0.023 * Re_v^(4/5) * Pr_v^0.4
    @test Nu_DittusBoelter(Re_v, Pr_v) == expected
  end

  @testset "f_Pethukov" begin
    Re_v = 10000.0; Pr_v = 0.7
    expected = (0.79 * log(Re_v) - 1.64)^(-2)
    @test f_Pethukov(Re_v, Pr_v) == expected
  end

  @testset "get_h_from_Nu" begin
    Nu = 100.0; k = 0.5; D = 0.1
    expected = Nu * k / D
    @test get_h_from_Nu(Nu, k, D) == expected
  end

  @testset "f_Haaland" begin
    Re_v = 10000.0; e_D = 0.001
    expected = (-1.8 * log10((e_D / 3.7)^1.11 + 6.9 / Re_v))^(-2)
    @test f_Haaland(Re_v, e_D) == expected
  end

  @testset "Schmidt" begin
    D = 0.1; mu = 0.001; rho = 1000.0
    expected = mu / (rho * D)
    @test Schmidt(D, mu, rho) == expected
  end

  @testset "Sherwood" begin
    Sc = 0.7; Re_v = 10000.0
    expected = 0.0096 * Re_v^0.913 * Sc^0.346
    @test Sherwood(Sc, Re_v) == expected
  end

  @testset "Sherwood_HT_analogy" begin
    Re_v = 10000.0; Sc = 0.7
    expected = 0.023 * Re_v^0.8 * Sc^0.4
    @test Sherwood_HT_analogy(Re_v, Sc) == expected
  end

  @testset "get_k_from_Sh" begin
    Sh = 100.0; L = 1.0; D = 0.1
    expected = Sh * D / L
    @test get_k_from_Sh(Sh, L, D) == expected
  end

  @testset "Re" begin
    rho = 1000.0; u = 2.0; L = 0.5; mu = 0.001
    expected = rho * u * L / mu
    @test Re(rho, u, L, mu) == expected
  end

  @testset "Pr" begin
    c_p = 1000.0; mu = 0.001; k = 0.5
    expected = mu * c_p / k
    @test Pr(c_p, mu, k) == expected
  end

  @testset "get_length_HX" begin
    deltaTML = 20.0; d_hyd = 0.1; U = 100.0; Q = 5000.0
    expected = Q / (U * π * d_hyd * deltaTML)
    @test get_length_HX(deltaTML, d_hyd, U, Q) == expected
  end

  @testset "Sherwood_bubbles" begin
    Sc = 0.7; Re_v = 10000.0
    expected = 0.089 * Re_v^0.69 * Sc^0.33
    @test Sherwood_bubbles(Sc, Re_v) == expected
  end

end
