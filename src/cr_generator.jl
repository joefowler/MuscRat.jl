using Statistics, Unitful

"""
    CRGenerator(....)

A cosmic ray particle generator.
"""
struct CRGenerator{P<:Unitful.Momentum, CRF<:CRFlux}
    Np::Int
    Nang::Int
    Pmin::P
    Pmax::P
    flux::CRF  # units are counts cm^-2 s^-1 (energy and solid angle already integrated over)

    logPGeVlim::LinRange
    cosθlim::LinRange
    relativeprob::Matrix{Float64}
    boxCDF::Vector{Float64}
    pij::Matrix{Float64}
end

function CRGenerator(
    Np::Integer, Nang::Integer, Pmin::P, Pmax::P,
    logPGeVlim::LinRange, cosθlim::LinRange, 
    Exspectrum::AbstractMatrix) where P<:Unitful.Momentum

    # Integrate s over all boxes
    relativeprob = zeros(eltype(Exspectrum), Np, Nang)
    for i=1:Np, j=1:Nang
        relativeprob[i,j] = mean(Exspectrum[i:i+1, j:j+1])
    end
    Δcos = cosθlim[2]-cosθlim[1]
    ΔlogP = logPGeVlim[2]-logPGeVlim[1]
    total_flux = 2π*Unitful.sr*sum(relativeprob)*(Δcos*ΔlogP)
    prob = relativeprob ./ sum(relativeprob)
    boxCDF = cumsum(vec(prob))  # runs through first column, then 2nd, etc.

    # Now create the 4 points for each box
    pij = Array{Float64}(undef, Np*Nang, 4)
    for j=1:Nang, i=1:Np
        k=Np*(j-1)+i
        rp = relativeprob[i,j]
        pij[k,1] = Exspectrum[i,j]/rp
        pij[k,2] = Exspectrum[i+1,j]/rp
        pij[k,3] = Exspectrum[i,j+1]/rp
        pij[k,4] = Exspectrum[i+1,j+1]/rp
    end
    CRGenerator(Np, Nang, Pmin, Pmax, total_flux, logPGeVlim, cosθlim, prob, boxCDF, pij)
end


"""
CRMuonGenerator(Np, Nang; Pmin=0.1, Pmax=1000, useReyna=false)

Returns a cosmic ray muon generator. The numerical integration is perfomed using `Np` momentum bins and
`Nang` angle bins. The momentum bins are logarithmically spaced between the energy `Pmin` and `Pmax`
(units: GeV/c). The angle bins are linearly spaced in cosθ between 0 and 1. The `useReyna` indicates whether
to use the Reyna 2006 spectrum model (the default is Charzidakis et al. 2015).

Use the generator to generate `N` muons like this:

```
generator = CRMuonGenerator(100, 100);
N = 1000000;
p,cosθ = generate(generator, N);
```
"""
function CRMuonGenerator(Np::Integer, Nang::Integer; 
    Pmin::P=0.1u"GeV/c", Pmax::P=1000.0u"GeV/c", useReyna::Bool=false
    ) where P<:Unitful.Momentum

    logPGeVlim = LinRange(log(Pmin/GeVc), log(Pmax/GeVc), 1+Np)
    cosθlim = LinRange(0, 1, 1+Nang)
    if useReyna
        spectrum = µspectrum_reyna_p
    else
        spectrum = µspectrum_chatzidakis_p
    end
    esval = spectrum(GeVc, 1)*Unitful.GeV
    s = zeros(eltype(esval), 1+Np, 1+Nang)
    for i=1:Np+1
        p = exp(logPGeVlim[i])*GeVc
        KE = muon_E(p)
        for j=1:Nang+1
            c = cosθlim[j]
            s[i,j] = KE*spectrum(p, c)
        end
    end
    CRGenerator(Np, Nang, Pmin, Pmax, logPGeVlim, cosθlim, s)
end



function CRGenerator(filename::AbstractString, mass::Real)
    lines = readlines(filename)
    header = split(lines[1], ",")
    Nang = parse(Int, header[1])
    Np = parse(Int, header[2])
    Pmin = parse(Float64, header[3])
    Pmax = parse(Float64, header[4])

    @assert Np == length(lines)-2
    @assert Nang == length(split(lines[2], ","))-2
    spectrum = zeros(Float64, 1+Np, 1+Nang)
    cosθlim = LinRange(0, 1, 1+Nang)
    ke = Float64[]
    for i = 1:Np+1
        energystr, spectrumstr... = split(lines[i+1], ",")
        push!(ke, parse(Float64, energystr))
        spectrum[i,:] .= [parse(Float64, x) for x in spectrumstr]
    end
    p = sqrt.((ke.+mass).^2 .- mass^2)
    dlogp = diff(log.(p))
    @show maximum(dlogp)-minimum(dlogp), mean(dlogp)*1e-3
    @assert maximum(dlogp)-minimum(dlogp) < mean(dlogp)*1e-3
    @assert abs(log(p[1]/Pmin)) < 1e-4
    @assert abs(log(p[end]/Pmax)) < 1e-4
    logPGeVlim = LinRange(log(Pmin), log(Pmax), length(p))
    for i=1:Nang+1
        spectrum[:,i] .*= p
    end

    CRGenerator(Np, Nang, Pmin, Pmax, logPGeVlim, cosθlim, spectrum)
end

@enum Particle begin
    Gamma
    µplus
    µminus
    Electron
    Positron
end


function ParmaGenerator(p::Particle)
    localpaths = Dict(
        Gamma => "data/parma_gamma.txt",
        µplus => "data/parma_mu+.txt",
        µminus => "data/parma_mu-.txt",
        Electron => "data/parma_electron.txt",
        Positron => "data/parma_positron.txt",
    )
    masses = Dict(
        Gamma => 0.0,
        µplus => mµ,
        µminus => mµ,
        Electron => me,
        Positron => me,
    )
    project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
    return CRGenerator(project_path(localpaths[p]), masses[p]*1000)
end


"""
    generate(mg::CRGenerator, N)

Return `(p,cosθ)`, two vectors of length `N`, that are randomly generated cosmic ray muons
using the generator `mg`. The momenta `p` are in units of GeV/c. The `cosθ` is the cosine of the zenith
angle (thus cosθ=1 means vertical, and 0 means horizontal).
"""
function generate(crg::CRGenerator, N::Integer)
    boxid = [findfirst(x->x>r, crg.boxCDF) for r in rand(N)]

    cosθ = Array{Float64}(undef, N)
    logp = Array{Float64}(undef, N)
    dθ = crg.cosθlim[2]-crg.cosθlim[1]
    dlogp = crg.logPGeVlim[2]-crg.logPGeVlim[1]
    for i=1:N
        k = boxid[i]
        boxi = 1 + ((k-1) % crg.Np)
        boxj = 1 + div(k-1, crg.Np)
        B = 0.5*(crg.pij[k,1] + crg.pij[k,2])
        A = 1 - B
        Fy = rand()
        y = (sqrt(B^2 + 4A*Fy)-B)/2A
        cosθ[i] = y*dθ + crg.cosθlim[boxj]

        Bx = (1-y)*crg.pij[k,1] + y*crg.pij[k,3]
        Ax = 0.5*((1-y)*(crg.pij[k,2]-crg.pij[k,1]) + y*(crg.pij[k,4]-crg.pij[k,3]))
        rescale = Ax+Bx  # replace conditional relativeprob with given y.
        Ax /= rescale
        Bx /= rescale
        Fx = rand()
        x = (sqrt(Bx^2 + 4Ax*Fx)-Bx)/2Ax
        logp[i] = x*dlogp + crg.logPGeVlim[boxi]
    end
    p = exp.(logp)*1u"GeV/c"
    p, cosθ
end

"""
    fluence(crg::CRGenerator, N=1)

Return the total cross-sectional area x time product that corresponds to `N` cosmic rays generated by
generator `crg`.
"""
fluence(crg::CRGenerator, N::Real=1) = N/crg.flux


function path_values(obj::Solid, cosθ::AbstractVector)
    tube_radius = smallest_radius(obj)
    N = length(cosθ)
    L = Array{Float64}(undef, N)
    for i=1:N
        # Random azimuthal angle
        ϕ = 2π*rand()
        sϕ = sin(ϕ)
        cϕ = cos(ϕ)
        cθ = cosθ[i]
        sθ = sqrt(1-cθ^2)
        n = [sθ*cϕ, sθ*sϕ, -cθ]

        # Random 2d position in the tube plane: polar coords (r,ψ) and cartesian (x0,y0)
        r = sqrt(rand())*tube_radius
        ψ = 2π*rand()
        x0 = r.*cos.(ψ)
        y0 = r.*sin.(ψ)

        # Where does the path cross the plane ⟂ the tube passing through the origin?
        # Rotate (x0,y0) by θ about the +x axis, then by ϕ about the z axis
        # x1, y1, z1 = x0, cθ*y0, sθ*y0
        # x2, y2, z2 = cϕ*x1+sϕ*y1, -sϕ*x1+cϕ*y1, z1
        origin = [cϕ*x0+sϕ*cθ*y0, -sϕ*x0+cθ*cϕ*y0, sθ*x0]
        # origin = [cϕ*cθ*x0-sϕ*y0, sϕ*cθ*x0+cϕ*y0, -sθ*x0]

        line = Line(n, origin)
        L[i] = path(obj, line)
    end
    L
end
