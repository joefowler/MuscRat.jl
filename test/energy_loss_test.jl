import Unitful.MeV
MeV_over_c = 1u"MeV/c"
MeV_per_cm = 1u"MeV/cm"

@testset "energy loss" begin
    tests = (
        (:Silicon, :µ,  300MeV_over_c, 3.9494MeV_per_cm),
        (:Copper,  :µ,  300MeV_over_c, 13.06MeV_per_cm),
        (:NaI,     :µ,  300MeV_over_c, 5.0MeV_per_cm),
        (:Silicon, :π,  300MeV_over_c, 4.06446MeV_per_cm),
        (:Silicon, :p, 3000MeV_over_c, 3.94683MeV_per_cm),
        (:Silicon, :e,    1MeV_over_c, 3.7292MeV_per_cm),
    )
    for t in tests
        material, particle, momentum, loss = t
        f = Eloss_function(material, particle)
        @test isapprox(f(momentum), loss; rtol=.001)
    end
    f = Eloss_function()
    @test f(300MeV_over_c) < f(1000MeV_over_c)
    @test f(300MeV_over_c) < f(100MeV_over_c)
end
