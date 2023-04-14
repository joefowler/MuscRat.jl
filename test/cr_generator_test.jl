@testset "muon generator" begin
    g = CRMuonGenerator(100, 100)
    N = 3000
    p,cosθ = generate(g, N)
    @test minimum(p) ≥ 0.1u"GeV/c"
    @test maximum(p) ≤ 1000u"GeV/c"
    @test minimum(cosθ) ≥ 0
    @test maximum(cosθ) ≤ 1

    Pmin = 1u"GeV/c"
    Pmax = 9u"GeV/c"
    g = CRMuonGenerator(100, 100; Pmin=Pmin, Pmax=Pmax)
    p,cosθ = generate(g, N)
    @test minimum(p) ≥ Pmin
    @test maximum(p) ≤ Pmax
    @test minimum(cosθ) ≥ 0
    @test maximum(cosθ) ≤ 1

    g1 = ParmaGenerator(MuscRat.µplus)
    g2 = ParmaGenerator(MuscRat.µplus; Pmin=2.0u"MeV/c")
    @test g1.Np > g2.Np
    @test_throws ErrorException g = ParmaGenerator(MuscRat.µplus; Pmax=1.0u"GeV/c", Pmin=10.0u"GeV/c")
    # Check integer values of Pmin... doesn't work yet!!!
    # ParmaGenerator(MuscRat.µplus; Pmin=2u"MeV/c")
end
