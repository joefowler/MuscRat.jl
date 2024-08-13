
@enum ParmaParticle begin
    ppNeutron = 0
    ppProton = 1
    ppHelium = 2
    ppµplus = 29
    ppµminus = 30
    ppElectron = 31
    ppPositron = 32
    ppGamma = 33
end

library = "parma/parma"

"""
    Parma_getHP(year::Integer, month::Integer, day::Integer)

Return the W-index of solar activity at date `year`-`month`-`day`.
"""
function Parma_getHP(year::Integer, month::Integer, day::Integer)
    Windex = @ccall library._Z8getHPcppiii(year::Cint, month::Cint, day::Cint)::Cdouble
end


"""
    Parma_getr(latitude::Real, longitude::Real=100)

Return the vertical cut-off rigidity (GV) at `latitude` and `longitude` (both in degrees).
"""
function Parma_getr(latitude::Real, longitude::Real=100)
    rigidity = @ccall library._Z7getrcppdd(latitude::Cdouble, longitude::Cdouble)::Cdouble
end

"""
    Parma_getd(altitude::Real, latitude::Real=100)

Return the atmospheric depth at `alititude` (in km above sea level) and at 
`latitude` (in degrees), in g/cm^2. 

If |`latitude`| ≤ 90, use NRLMSISE-00 data.
Otherwise, use the U.S. Standard Atmosphere.
"""
function Parma_getd(altitude::Real, latitude::Real=100)
    depth = @ccall library._Z7getdcppdd(altitude::Cdouble, latitude::Cdouble)::Cdouble
end

"""
    Parma_getSpec(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real)

Get the cosmic ray spectrum for particle `particleID` at `energy` (in MeV per nucleon), also given the
W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV), the atmospheric depth `atmosphericDepth` (g/cm^2)
and the local geometry parameter.

`geometryParam` should be in [0, 1] if a water weight-fraction, or 10 for no earth, or 100 for the "black hole" model,
or in [-10, 0] for pilot, and g<-10 for cabin. I don't know what that all means. Use 0.15.
"""
function Parma_getSpec(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real)
    spectrum = @ccall library._Z10getSpecCppiddddd(particleID::Cint, Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble, energy::Cdouble, geometryParam::Cdouble)::Cdouble
end

"""
    Parma_getSpecAngFinal(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::Real)

Get the angular cosmic ray spectrum for particle `particleID` at `energy` (in MeV per nucleon) and zenith
angle `cosθ`, also given the W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV), the atmospheric
depth `atmosphericDepth` (g/cm^2) and the local geometry parameter.

`geometryParam` should be in [0, 1] if a water weight-fraction, or 10 for no earth, or 100 for the "black hole" model,
or in [-10, 0] for pilot, and g<-10 for cabin. I don't know what that all means. Use 0.15.
"""
function Parma_getSpecAngFinal(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::Real)
    iangParticle = Dict(ppNeutron=>1, ppProton=>2, ppHelium=>3, ppµminus=>4, ppµplus=>4, ppElectron=>5, ppPositron=>5, ppGamma=>6)
    id = get(iangParticle, particleID, 0)
    spectrum = @ccall library._Z18getSpecAngFinalCppidddddd(id::Cint, Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble, energy::Cdouble, geometryParam::Cdouble, cosθ::Cdouble)::Cdouble
end

"""
    Parma_get511flux(Windex::Real, cutoffRigidity::Real, atmosphericDepth::Real)

Get the flux of 511 keV photons (cm^-2 s^-1), given the W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV),
and the atmospheric depth `atmosphericDepth` (g/cm^2).
"""
function Parma_get511flux(Windex::Real, cutoffRigidity::Real, atmosphericDepth::Real)
    spectrum = @ccall library._Z13get511fluxCppddd(Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble)::Cdouble
end

