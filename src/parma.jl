using FastGaussQuadrature

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

function parmaID(p::Particle)
    p == Gamma && return ppGamma
    p == µplus && return ppµplus
    p == µminus && return ppµminus
    p == Electron && return ppElectron
    p == Positron && return ppPositron
    p == Neutron && return ppNeutron
    p == Proton && return ppProton
    nothing
end

const ParmaCompiledLibrary = "parma/parma"

"""
    Parma_getHP(year::Integer, month::Integer, day::Integer)

Return the W-index of solar activity at date `year`-`month`-`day`.
"""
function Parma_getHP(year::Integer, month::Integer, day::Integer)
    Windex = @ccall ParmaCompiledLibrary._Z8getHPcppiii(year::Cint, month::Cint, day::Cint)::Cdouble
end


"""
    Parma_getr(latitude::Real, longitude::Real=100)

Return the vertical cut-off rigidity (GV) at `latitude` and `longitude` (both in degrees).
"""
function Parma_getr(latitude::Real, longitude::Real=100)
    rigidity = @ccall ParmaCompiledLibrary._Z7getrcppdd(latitude::Cdouble, longitude::Cdouble)::Cdouble
end

"""
    Parma_depth(altitude::Real, latitude::Real=100)

Return the atmospheric depth at `alititude` (in km above sea level) and at 
`latitude` (in degrees), in g/cm^2. 

If |`latitude`| ≤ 90, use NRLMSISE-00 data.
Otherwise, use the U.S. Standard Atmosphere.
"""
function Parma_depth(altitude::Real, latitude::Real=100)
    depth = @ccall ParmaCompiledLibrary._Z7getdcppdd(altitude::Cdouble, latitude::Cdouble)::Cdouble
end

"""
    Parma_getSpec(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real)

Get the angle-integrated cosmic ray spectrum for particle `particleID` at `energy` (in MeV per nucleon), also given the
W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV), the atmospheric depth `atmosphericDepth` (g/cm^2)
and the local geometry parameter `geometryParam`. 

The angle-integrated spectrum's units are events per (cm^2⋅s⋅MeV/n).

`geometryParam` should be in [0, 1] if a water weight-fraction of nearby earth, or 10 for no earth, or 100 for the "black hole" model,
or in [-10, 0] for pilot, and g<-10 for cabin. I don't know what that all means. Use 0.15.
"""
function Parma_getSpec(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real)
    spectrum = @ccall ParmaCompiledLibrary._Z10getSpecCppiddddd(particleID::Cint, Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble, energy::Cdouble, geometryParam::Cdouble)::Cdouble
end

"""
    Parma_getSpecAngFinal(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::Real)

Get the angular differential cosmic ray spectrum for particle `particleID` at `energy` (in MeV per nucleon) and zenith
angle `cosθ`, also given the W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV), the atmospheric
depth `atmosphericDepth` (g/cm^2) and the local geometry parameter.

The angular differential spectrum units are sr^{-1}. This is a factor meant to multiply the "omnidirectional" spectrum.

WARNING: It is designed to have an integral over all 4π sr equal to 1, but owing to the details of the smoothing over various
parameters over energy, it will not actually accomplish this. It's on you to renormalize it appropriately.

`geometryParam` should be in [0, 1] if a water weight-fraction, or 10 for no earth, or 100 for the "black hole" model,
or in [-10, 0] for pilot, and g<-10 for cabin. I don't know what that all means. Use 0.15.
"""
function Parma_getSpecAngFinal(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::Real)
    iangParticle = Dict(ppNeutron=>1, ppProton=>2, ppHelium=>3, ppµminus=>4, ppµplus=>4, ppElectron=>5, ppPositron=>5, ppGamma=>6)
    id = get(iangParticle, particleID, 0)
    angular_weighting = @ccall ParmaCompiledLibrary._Z18getSpecAngFinalCppidddddd(id::Cint, Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble, energy::Cdouble, geometryParam::Cdouble, cosθ::Cdouble)::Cdouble
end

function Parma_angdist(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::Real)
    f(x) = Parma_getSpecAngFinal(particleID, Windex, cutoffRigidity, atmosphericDepth, energy, geometryParam, x)
    
    x, w = gausslegendre(150)
    integral = 2π * dot(f.(x), w)
    f(cosθ)/integral
end

function Parma_angdist(particleID::ParmaParticle, Windex::Real, cutoffRigidity::Real, 
    atmosphericDepth::Real, energy::Real, geometryParam::Real, cosθ::AbstractArray)
    f(x) = Parma_getSpecAngFinal(particleID, Windex, cutoffRigidity, atmosphericDepth, energy, geometryParam, x)
    
    x, w = gausslegendre(150)
    integral = 2π * dot(f.(x), w)
    f.(cosθ)/integral
end

"""
    Parma_get511flux(Windex::Real, cutoffRigidity::Real, atmosphericDepth::Real)

Get the flux of 511 keV photons (cm^-2 s^-1), given the W-index `Windex`, the cutoff rigidity `cutoffRigidity` (in GV),
and the atmospheric depth `atmosphericDepth` (g/cm^2).
"""
function Parma_get511flux(Windex::Real, cutoffRigidity::Real, atmosphericDepth::Real)
    spectrum = @ccall ParmaCompiledLibrary._Z13get511fluxCppddd(Windex::Cdouble, cutoffRigidity::Cdouble, 
        atmosphericDepth::Cdouble)::Cdouble
end


struct CRObserver
    year::Int
    month::Int
    day::Int
    latitude::Float64
    longitude::Float64
    altitude::Float64
    g::Float64
    Windex::Float64
    cutoffRigidity::Float64
    depth::Float64
    function CRObserver(year,month,day,lat,long,alt,g=0.15)
        Windex = Parma_getHP(year, month, day)
        r = Parma_getr(lat, long)
        depth = Parma_depth(alt)
        new(year,month,day,lat,long,alt,g,Windex,r,depth)
    end
end

function CRspectrum(energy::Real, obs::CRObserver, id::Particle; depth::Real=obs.depth)
    pid = parmaID(id)
    Parma_getSpec(pid, obs.Windex, obs.cutoffRigidity, depth, energy, obs.g)
end

function CRspectrum(energy::AbstractArray, obs::CRObserver, id::Particle; depth::Real=obs.depth)
    pid = parmaID(id)
    f(e) = Parma_getSpec(pid, obs.Windex, obs.cutoffRigidity, depth, e, obs.g)
    f.(energy)
end

function CRangularSpectrum(energy::Real, cosθ::Real, obs::CRObserver, id::Particle; depth::Real=obs.depth)
    pid = parmaID(id)
    Parma_angdist(pid, obs.Windex, obs.cutoffRigidity, depth, energy, obs.g, cosθ) * CRspectrum(energy, obs, id, depth=depth)
end

function CRangularSpectrum(energy::Real, cosθ::AbstractArray, obs::CRObserver, id::Particle; depth::Real=obs.depth)
    pid = parmaID(id)
    Parma_angdist(pid, obs.Windex, obs.cutoffRigidity, depth, energy, obs.g, cosθ) * CRspectrum(energy, obs, id, depth=depth)
end

