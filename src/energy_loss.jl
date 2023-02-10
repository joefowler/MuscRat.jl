"""
    Eloss_function(material=:Silicon, particle=:µ)

Returns a _function_ `dEdX(p)` that accepts a momentum (MeV/c units) and returns the mean
energy loss according to the Bethe-Bloch formula (Eq 33.5 and 33.4 in the 2018 _Review of 
Particle Properties_).

* `material` is the absorbing material. Allowed values: (:Silicon, :Copper, :NaI)
* `particle` is the particle being stopped. Allowed values: (:µ, :e, :p)
"""
function Eloss_function(material=:Silicon, particle=:µ)
    masses = Dict(
        :µ => 105.65837,  # MeV/c^2
        :π => 139.5706,   # MeV/c^2
        :e => 0.5109989,  # MeV/c^2
        :p => 938.27208,  # MeV/c^2
    )
    M = masses[particle]
    K = 0.307075     # MeV g^-1 cm^2
    if material == :Silicon
        I = 175e-6  # silicon excitation energy (eV, not MeV, for some reason!) See Fig 33.5
        ZoverA = 14/28.085  # silicon charge-to-mass ratio
        ρ_absorber = 2.329 # g/cm^3
    elseif material == :Copper
        I = 29*11e-6
        ZoverA = 0.45636
        ρ_absorber = 8.960
    elseif material == :NaI
        I = .5*(140+500)*1e-6  # roughly 140 eV and 500 eV for Na and I separately...so average them?
        ZoverA = 0.42697
        ρ_absorber = 3.667
    else
        error("material must be one of (:Silicon, :Copper, :NaI)")
    end

    # This is the dE/dx value in units of MeV/cm
    function dEdx(p::Real)
        βγ = p/M
        γ = sqrt(1+βγ^2)
        β = βγ/γ
        mratio = masses[:e]/M
        Wmax_factor = 1/sqrt(1+2γ*mratio+mratio^2)
        bracket = log(2masses[:e]*βγ^2*Wmax_factor/I)-β^2
        ρ_absorber*K*(ZoverA)*bracket/β^2
    end
end
