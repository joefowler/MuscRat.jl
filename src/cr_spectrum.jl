# See README for full citations of Su, Charzidakis, and Reyna

function µspectrum_chatzidakis(Eµ, cosθ)
    ###############################################
    #
    #    Phenomenological model constants
    #
    A=0.002382;                # constant A
    λ=120;                     # absorption mean free path 120 g/cm2
    kappa=2.645;               # exponent (-)
    bp=0.771; jp=148.16;       # correction factor (-); factor (GeV)
    α = 0.0025;                # muon energy loss in GeV/g/cm2
    ρ = 0.76;                  # fraction of pion energy that is transferred to muon
    y0=1000;                   # atmoshperic depth g/cm2
    Bm=1.041231831;            # correction factor (-); 
    ###############################################
    Eµ < 0 && error("Negative µ kinetc energy")
    cosθ ≤ 0 && return 0.0
    secθ = 1/cosθ

    # Su eq (7): energy of the π that produced the µ
    ρEp1 = (Eµ+α*y0*(secθ-0.1))

    # Su eq (8): probability for µ to reach sea level
    Pm1 = (0.1*cosθ*(1 - α*(y0*secθ-100)/ρEp1))^(Bm/((ρEp1+100α)*cosθ))

    # Su eq (6): cosmic ray µ spectrum
    A*Pm1*((ρEp1/ρ)^(-kappa))*λ*bp*jp/(ρEp1/ρ*cosθ+bp*jp)
end

function µspectrum_reyna_p(pµ, cosθ)
    # The Bugaev 1998 model as revised by Reyna 2006 for non-vertical muons, and as
    # reported in Su et al. 2021.
    pµ < 0 && error("Negative µ momentum")
    cosθ ≤ 0 && return 0.0
    AB = 0.00253
    a = [0.2455, 1.288, -0.2555, 0.0209]
    y = log10(pµ*cosθ)   
    exponent = a[4]*y^3+a[3]*y^2+a[2]*y+a[1]
    ΦB =  AB * (pµ*cosθ)^(-exponent)
    cosθ^3 * ΦB
end

const mµ = 0.105659  # µ mass in GeV
const me = 0.000510999 # e± mass in GeV

muon_E(p) = p<0 ? error("Negative µ momentum") : sqrt(p^2+mµ^2)-mµ
muon_p(E) = E<0 ? error("Negative µ kinetc energy") : sqrt((E+mµ)^2-mµ^2)

µspectrum_chatzidakis_p(p, cosθ) = µspectrum_chatzidakis(muon_E(p), cosθ)
µspectrum_reyna(E, cosθ) = µspectrum_reyna_p(muon_p(E), cosθ)

"""
    µspectrum(E, cosθ)

Return the differential cosmic ray muon energy spectrum in [counts GeV^-1 s^-1 sr^-1 cm^-2].
- `E`    is the µ kinetic energy in GeV
- `cosθ` is the cosine of the zenith angle
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
