using Unitful

# See README for full citations of Su, Charzidakis, and Reyna
function µspectrum_chatzidakis(Eµ, cosθ; y=1000u"g/cm^2")
    ###############################################
    #
    #    Phenomenological model constants
    #
    # y = atmospheric depth g/cm2
    A = 0.002382u"1/GeV/g/s/sr"; # constant A, a fitting parameter of model
    λ = 120u"g/cm^2";            # absorption mean free path 120 g/cm2
    kappa = 2.645;               # exponent, a fitting parameter of model
    b = 0.771;                   # correction factor (improves on an isothermal atmosphere assumption)
    jπ = 148.16u"GeV";           # π j-factor(GeV)
    α = 2.5u"MeV/(g/cm^2)";      # µ energy loss in MeV/g/cm2
    ρ = 0.76;                    # fraction of π energy that is transferred to µ
    Bµ = 1.041231831u"GeV";      # correction factor (-);
    y0 = 1000u"g/cm^2"           # atmospheric depth at sea level
    Bµ *= y/y0
    ###############################################
    Eµ < Ezero && error("Negative µ kinetic energy")
    cosθ ≤ 0 && return 0.0A*λ
    secθ = 1/cosθ

    # Su eq (7): energy of the π that produced the µ
    ρEπ = (Eµ+α*y*(secθ-0.1))
    Eπ = ρEπ/ρ

    # Su eq (8): probability for µ to reach sea level
    scale100 = 100u"g/cm^2"
    power8 = convert(Float64, Bµ/((ρEπ+scale100*α)*cosθ))
    Pµ1 = (0.1*cosθ*(y/y0)*(1 - α*(y*secθ-scale100)/ρEπ))^power8

    # Su eq (6): cosmic ray µ spectrum
    A*Pµ1*((Eπ/1u"GeV" |> NoUnits)^(-kappa))*λ*(b*jπ/(Eπ*cosθ+b*jπ) |> NoUnits)
end

function µspectrum_reyna_p(pµ, cosθ)
    # The Bugaev 1998 model as revised by Reyna 2006 for non-vertical muons, and as
    # reported in Su et al. 2021.
    pµ < Pzero && error("Negative µ momentum")
    AB = 0.00253*u"1/GeV/s/sr/cm^2"
    cosθ ≤ 0 && return 0AB
    a = [0.2455, 1.288, -0.2555, 0.0209]
    pµGeV = convert(Float64, pµ/u"GeV/c")
    y = log10(pµGeV*cosθ)
    exponent = a[4]*y^3+a[3]*y^2+a[2]*y+a[1]
    ΦBvert =  AB * (pµGeV*cosθ)^(-exponent)
    cosθ^3 * ΦBvert
end

muon_E(p) = p < Pzero ? error("Negative µ momentum") : sqrt((p*Unitful.c)^2+(mµ*Unitful.c^2)^2)-mµ*Unitful.c^2
muon_p(E) = E < Ezero ? error("Negative µ kinetic energy") : sqrt((E+mµ*Unitful.c^2)^2-mµ^2*Unitful.c^4)/Unitful.c

µspectrum_chatzidakis_p(p, cosθ; y=1000u"g/cm^2") = µspectrum_chatzidakis(muon_E(p), cosθ; y=y)
µspectrum_reyna(E, cosθ) = µspectrum_reyna_p(muon_p(E), cosθ)

"""
    µspectrum(E, cosθ; y0=1000u"g/cm^2")

Return the differential cosmic ray muon energy spectrum in [counts GeV^-1 s^-1 sr^-1 cm^-2].
- `E`    is the µ kinetic energy in GeV
- `cosθ` is the cosine of the zenith angle
- `y0`   is the atmospheric overburden (sea level is approximately 1000 g/cm^2)
Uses the Chatzidakis 2015 spectrum (doi:10.1016/j.nima.2015.09.033)
See `µspectrum_reyna` to use the Reyna 2006 spectrum (arXiv:hep-ph/0604145)
See `µspectrum_p` or `µspectrum_reyna_p` to accept momentum instead of kinetic energy as the 1st parameter
"""
µspectrum = µspectrum_chatzidakis

"""
    µspectrum_p(p, cosθ)

Return the differential cosmic ray muon energy spectrum in [counts GeV^-1 s^-1 sr^-1 cm^-2].
- `p`    is the µ momentum in GeV/c
- `cosθ` is the cosine of the zenith angle
Uses the Chatzidakis 2015 spectrum (doi:10.1016/j.nima.2015.09.033)
See `µspectrum_reyna_p` to use the Reyna 2006 spectrum (arXiv:hep-ph/0604145)
See `µspectrum` or `µspectrum_reyna` to accept kinetic energy instead of momentum as the 1st parameter
"""
µspectrum_p = µspectrum_chatzidakis_p
