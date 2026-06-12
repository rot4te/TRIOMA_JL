@testset "Hydraulic diameter — cross-section types" begin

    @testset "Circular" begin
        cs = Circular(0.3)
        @test flow_area(cs) ≈ π * 0.15^2
        @test wetted_perimeter(cs) ≈ π * 0.3
        @test hydraulic_diameter(cs) ≈ 0.3
    end

    @testset "Rectangular — square" begin
        cs = Rectangular(0.2, 0.2)
        @test flow_area(cs) ≈ 0.04
        @test wetted_perimeter(cs) ≈ 0.8
        @test hydraulic_diameter(cs) ≈ 0.2
    end

    @testset "Rectangular — 2:1 aspect ratio" begin
        cs = Rectangular(0.2, 0.1)
        @test flow_area(cs) ≈ 0.02
        @test wetted_perimeter(cs) ≈ 0.6
        @test hydraulic_diameter(cs) ≈ 4 * 0.02 / 0.6
    end

    @testset "Annulus" begin
        cs = Annulus(0.3, 0.1)
        @test flow_area(cs) ≈ π * (0.15^2 - 0.05^2)
        @test wetted_perimeter(cs) ≈ π * (0.3 + 0.1)
        # D_h = (D_outer^2 - D_inner^2) / (D_outer + D_inner) = D_outer - D_inner
        @test hydraulic_diameter(cs) ≈ 0.2
    end

    @testset "TwistedElliptical — circular degenerate (a = b)" begin
        cs = TwistedElliptical(0.15, 0.15)
        @test flow_area(cs) ≈ π * 0.15^2
        @test wetted_perimeter(cs) ≈ π * 0.3
        @test hydraulic_diameter(cs) ≈ 0.3
    end

    @testset "TwistedElliptical — 2:1 aspect ratio" begin
        cs = TwistedElliptical(0.2, 0.1)
        @test flow_area(cs) ≈ π * 0.02
        @test hydraulic_diameter(cs) ≈ 4 * flow_area(cs) / wetted_perimeter(cs)
        @test hydraulic_diameter(cs) > 0.1
        @test hydraulic_diameter(cs) < 0.4
    end

    @testset "TwistedElliptical — high aspect ratio (5:1)" begin
        cs = TwistedElliptical(0.25, 0.05)
        @test flow_area(cs) ≈ π * 0.0125
        @test hydraulic_diameter(cs) ≈ 4 * flow_area(cs) / wetted_perimeter(cs)
    end

    @testset "Geometry — circular fallback" begin
        g = Geometry(L=1.0, D=0.3, dw=0.01)
        @test hydraulic_diameter(g) ≈ 0.3
        @test flow_area(g) ≈ π * 0.15^2
        @test wetted_perimeter(g) ≈ π * 0.3
    end

    @testset "Geometry — explicit TwistedElliptical cross-section" begin
        cs = TwistedElliptical(0.2, 0.1)
        g  = Geometry(L=1.0, D=0.3, dw=0.01, cross_section=cs)
        @test flow_area(g)          ≈ flow_area(cs)
        @test wetted_perimeter(g)   ≈ wetted_perimeter(cs)
        @test hydraulic_diameter(g) ≈ hydraulic_diameter(cs)
        @test get_fluid_volume(g)   ≈ flow_area(cs) * 1.0
    end

    @testset "define_component_volumes! auto-derives d_Hyd for TwistedElliptical" begin
        cs   = TwistedElliptical(0.2, 0.1)
        comp = Component(
            geometry = Geometry(L=1.0, D=0.3, dw=0.01, cross_section=cs),
            fluid    = Fluid(MS=true),
            membrane = Membrane()
        )
        define_component_volumes!(comp)
        @test comp.fluid.V     ≈ flow_area(cs) * 1.0
        @test comp.fluid.d_Hyd ≈ hydraulic_diameter(cs)
    end

end

@testset "Twisted-tube correlations (tube-side)" begin

    Re_t = 50_000.0; Pr_t = 7.0
    D_max = 0.022; D_min = 0.009; s = 0.2
    D_h   = hydraulic_diameter(TwistedElliptical(D_max/2, D_min/2))

    @testset "Nu_Iev_tube" begin
        Nu = Nu_Iev_tube(Re_t, D_max, s)
        @test Nu ≈ 0.019 * Re_t^0.8 * (1 + 0.547 * (s/D_max)^(-0.83))
        # Tighter twist (smaller s) → higher Nu
        @test Nu_Iev_tube(Re_t, D_max, s/2) > Nu
    end

    @testset "Nu_Asma_tube" begin
        Nu = Nu_Asma_tube(Re_t, Pr_t, D_max, s)
        @test Nu ≈ 0.021 * Re_t^0.8 * Pr_t^0.4 * (1 + 2.1 * (s/D_max)^(-0.91))
        @test Nu_Asma_tube(Re_t, Pr_t, D_max, s/2) > Nu
    end

    @testset "Nu_Si_tube" begin
        Nu = Nu_Si_tube(Re_t, Pr_t, D_h, D_max, s)
        @test Nu ≈ 0.396 * Re_t^0.544 * (s/D_h)^0.161 * (s/D_max)^(-0.519) * Pr_t^0.33
    end

    @testset "Nu_Yang1_tube" begin
        Nu = Nu_Yang1_tube(Re_t, Pr_t, D_h, D_max, D_min, s)
        @test Nu ≈ 0.034 * Re_t^0.784 * Pr_t^0.333 *
                   (D_min/D_max)^(-0.590) * (s/D_h)^(-0.165)
    end

    @testset "Nu_Yang2_tube" begin
        Nu = Nu_Yang2_tube(Re_t, Pr_t, D_h, D_max, D_min, s)
        @test Nu ≈ 1.50618 * Re_t^0.51825 * Pr_t^(-1.2446) *
                   (D_max/D_min)^1.12252 * (s/D_h)^(-0.32367)
    end

    @testset "fD_Iev_tube" begin
        f = fD_Iev_tube(Re_t, D_max, s)
        @test f ≈ 0.316 * (1 + 3.27 * (s/D_max)^(-0.87)) * Re_t^(-0.25)
        # Tighter twist → higher friction
        @test fD_Iev_tube(Re_t, D_max, s/2) > f
    end

    @testset "fD_Asma_tube" begin
        f = fD_Asma_tube(Re_t, D_max, s)
        @test f ≈ 0.82 * (s/D_max)^(-0.63) * Re_t^(-0.18)
        @test fD_Asma_tube(Re_t, D_max, s/2) > f
    end

    @testset "fD_Si_tube" begin
        f = fD_Si_tube(Re_t, D_max, s)
        r  = s / D_max
        a1 = -19.70 + 4.90*r - 0.22*r^2
        a2 =  10.52 - 2.66*r + 0.12*r^2
        a3 =  -1.47 + 0.36*r - 0.016*r^2
        # log10, confirmed by Hughes (2017) Appendix C Python code
        @test f ≈ 10^(a1 + a2*log10(Re_t) + a3*log10(Re_t)^2)
        @test f > 0
    end

    @testset "fD_Gao_tube" begin
        f = fD_Gao_tube(Re_t, D_h, D_max, D_min, s)
        @test f ≈ 4.572 * Re_t^(-0.521) * (D_min/D_max)^(-0.334) * (s/D_h)^(-0.082)
        @test f > 0
    end

    @testset "fD_Yang_tube" begin
        f = fD_Yang_tube(Re_t, Pr_t, D_h, D_max, D_min, s)
        @test f ≈ 0.71497 * Re_t^0.07777 * Pr_t^(-1.03974) *
                  (D_max/D_min)^(-0.076212) * (s/D_h)^(-0.33393)
        @test f > 0
    end

    @testset "Nu_Yang_lam" begin
        # Use laminar-range Re
        Re_lam = 200.0; Pr_lam = 7.0
        Nu = Nu_Yang_lam(Re_lam, Pr_lam, D_max, D_min, s)
        @test Nu ≈ 3.66 + 0.512 * Re_lam^0.477 * Pr_lam^0.975 *
                   (1 - D_min/D_max)^1.532 * (s/D_min)^(-0.609)
        # More eccentric tube (smaller D_min/D_max) → higher Nu
        @test Nu_Yang_lam(Re_lam, Pr_lam, D_max, D_min*0.5, s) > Nu
    end

end
