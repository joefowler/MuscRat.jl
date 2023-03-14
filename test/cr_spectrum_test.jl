using MuscRat, Unitful

@testset "spectrum" begin
    E = 1.0u"GeV"
    p = 1.0u"GeV/c"
    all_functions = (µspectrum, µspectrum_reyna, µspectrum_p, µspectrum_reyna_p)
    args = (E, E, p, p)
    for (f, arg) in zip(all_functions, args)
        @test ustrip(f(arg, -.3)) == 0.0
        @test_throws Exception f(-arg, .5)
    end
    mµ = MuscRat.mµ
    p = MuscRat.muon_p(E)
    @test µspectrum(E, 1) ≈ µspectrum_p(p, 1)
    @test µspectrum_reyna(E, 1) ≈ µspectrum_reyna_p(p, 1)
    @test isapprox(µspectrum(E, 1), µspectrum_reyna(E, 1); rtol = 0.30)
end
