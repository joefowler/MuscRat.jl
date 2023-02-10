@testset "muon generator" begin
    g = CRMuonGenerator(100, 100)
    N = 10000
    p,cosθ = generate(g, N)
    @test minimum(p) ≥ 0.1
    @test maximum(p) ≤ 1000
    @test minimum(cosθ) ≥ 0
    @test maximum(cosθ) ≤ 1
end
