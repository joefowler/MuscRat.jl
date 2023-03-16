import Unitful.MeV, Unitful.eV

"""
    Eloss_function(material=:Silicon, particle=:µ)

Returns a _function_ `dEdX(p)` that accepts a momentum (MeV/c units) and returns the mean
energy loss according to the Bethe-Bloch formula (Eq 33.5 and 33.4 in the 2018 _Review of
Particle Properties_).

* `material` is the absorbing material. Allowed values: (:Silicon, :Copper, :NaI)
* `particle` is the particle being stopped. Allowed values: (:µ, :e, :p)
"""
function Eloss_function(material=:Silicon, particle=:µ)
    massenergies = Dict(
        :µ => 105.65837MeV,
        :π => 139.5706MeV,
        :e => 0.5109989MeV,
        :p => 938.27208MeV,
    )
    Mc2 = massenergies[particle]
    K = 0.307075u"MeV/g*cm^2"
    if material == :Silicon
        I = 175eV  # silicon excitation energy (eV, not MeV, for some reason!) See Fig 33.5
        ZoverA = 14/28.085  # silicon charge-to-mass ratio
        ρ_absorber = 2.329u"g/cm^3"
    elseif material == :Copper
        I = 29*11eV
        ZoverA = 0.45636
        ρ_absorber = 8.960u"g/cm^3"
    elseif material == :NaI
        I = .5*(140eV+500eV)  # roughly 140 eV and 500 eV for Na and I separately...so average them?
        ZoverA = 0.42697
        ρ_absorber = 3.667u"g/cm^3"
    else
        error("material must be one of (:Silicon, :Copper, :NaI)")
    end

    # This is the dE/dx value in units of MeV/cm
    function dEdx(p::Unitful.Momentum)
        βγ = p*Unitful.c/Mc2
        γ = sqrt(1+βγ^2)
        β = βγ/γ
        mratio = massenergies[:e]/Mc2
        Wmax_factor = 1/sqrt(1+2γ*mratio+mratio^2)
        bracket = log(2massenergies[:e]*βγ^2*Wmax_factor/I)-β^2
        dedx = ρ_absorber*K*(ZoverA)*bracket/β^2
        uconvert(u"MeV/cm", dedx)
    end
end
