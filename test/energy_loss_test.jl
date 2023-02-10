using MuscRat

@testset "energy loss" begin
    tests = (
        (:Silicon, :µ, 300, 3.9494),
        (:Copper, :µ, 300, 13.06),
        (:NaI, :µ, 300, 5.0),
        (:Silicon, :π, 300, 4.06446),
        (:Silicon, :p, 3000, 3.94683),
        (:Silicon, :e, 1, 3.7292),
    )
    for t in tests
        material, particle, p_MeV, loss = t
        f = Eloss_function(material, particle)
        @test isapprox(f(p_MeV), loss; rtol=.001)
    end
    f = Eloss_function()
    @test f(300) < f(1000)
    @test f(300) < f(100)
end
