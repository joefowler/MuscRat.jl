import Unitful.MeV, Unitful.keV
MeV_over_c = 1u"MeV/c"
MeV_per_cm = 1u"MeV/cm"

@testset "energy loss" begin
    tests = (
        (:Silicon, :µ,  300MeV_over_c, 3.9494MeV_per_cm, 296.7keV),
        (:Copper,  :µ,  300MeV_over_c, 13.06MeV_per_cm, 1049.27keV),
        (:NaI,     :µ,  300MeV_over_c, 5.0MeV_per_cm, 375.66keV),
        (:Silicon, :π,  300MeV_over_c, 4.06446MeV_per_cm, 312.25keV),
        (:Silicon, :p, 3000MeV_over_c, 3.94683MeV_per_cm, 293.6keV),
        (:Silicon, :e,    1MeV_over_c, 3.7292MeV_per_cm, 320.97keV),
    )
    thick = 1u"mm"
    for t in tests
        material, particle, momentum, lossrate, loss = t
        (dEdx, Δprob) = Eloss_functions(material, particle)
        @test isapprox(dEdx(momentum), lossrate; rtol=.001)
        @test isapprox(Δprob(momentum, thick), loss; rtol=.001)
    end
    dEdx,Δprob = Eloss_functions()
    pminI = 350MeV_over_c
    @test dEdx(pminI) < dEdx(3pminI)
    @test dEdx(pminI) < dEdx(0.3pminI)
    @test Δprob(pminI, thick) < dEdx(pminI)*thick
    @test Δprob(pminI, thick) < Δprob(0.3pminI, thick)
end
