using MuscRat

@testset "spectrum" begin
    all_functions = (µspectrum, µspectrum_p, µspectrum_reyna, µspectrum_reyna_p)
    for f in all_functions
        @test f(1.0, -.3) == 0.0
        @test_throws Exception f(-1.0, .5) == 0.0
    end
    mµ = MuscRat.mµ
    E = 1.0
    p = sqrt((E+MuscRat.mµ)^2-mµ^2)
    @test µspectrum(E, 1) ≈ µspectrum_p(p, 1)
    @test µspectrum_reyna(E, 1) ≈ µspectrum_reyna_p(p, 1)
    @test isapprox(µspectrum(E, 1), µspectrum_reyna(E, 1); rtol = 0.30)
end
