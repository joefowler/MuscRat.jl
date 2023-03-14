@testset "muon generator" begin
    g = CRMuonGenerator(100, 100)
    N = 10000
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
end
