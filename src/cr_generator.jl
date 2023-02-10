using Statistics

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
struct CRMuonGenerator
    Np::Int
    Nang::Int
    Pmin::Float64
    Pmax::Float64
    flux::Float64  # units are counts cm^-2 s^-1 (energy and solid angle already integrated over)

    logPlim::LinRange
    cosθlim::LinRange
    prob::Matrix{Float64}
    boxCDF::Vector{Float64}
    pij::Matrix{Float64}

    function CRMuonGenerator(Np::Integer, Nang::Integer; Pmin::Real=0.1, Pmax::Real=1000.0, useReyna::Bool=false)
        logPlim = LinRange(log(Pmin), log(Pmax), 1+Np)
        cosθlim = LinRange(0, 1, 1+Nang)
        if useReyna
            spectrum = µspectrum_reyna_p
        else
            spectrum = µspectrum_chatzidakis_p
        end
        prob = zeros(Float64, Np, Nang)
        s = zeros(Float64, 1+Np, 1+Nang)
        for i=1:Np+1
            p = exp(logPlim[i])
            for j=1:Nang+1
                c = cosθlim[j]
                s[i,j] = p*spectrum(p, c)
            end
        end

        # Integrate s over all boxes
        for i=1:Np, j=1:Nang
            prob[i,j] = mean(s[i:i+1, j:j+1])
        end
        Δcos = cosθlim[2]-cosθlim[1]
        ΔlogP = logPlim[2]-logPlim[1]
        total_flux = 2π*sum(prob)*(Δcos*ΔlogP)
        prob ./= sum(prob)
        boxCDF = cumsum(vec(prob))  # runs through first column, then 2nd, etc.

        # Now create the 4 points for each box
        pij = Array{Float64}(undef, Np*Nang, 4)
        for j=1:Nang, i=1:Np
            k=Np*(j-1)+i
            pij[k,1] = s[i,j]
            pij[k,2] = s[i+1,j]
            pij[k,3] = s[i,j+1]
            pij[k,4] = s[i+1,j+1]
            pij[k,:] ./= mean(pij[k,:])
        end

        new(Np, Nang, Pmin, Pmax, total_flux, logPlim, cosθlim, prob, boxCDF, pij)
    end
end

"""
    generate(mg::CRMuonGenerator, N)

Return `(p,cosθ)`, two vectors of length `N`, that are randomly generated cosmic ray muons
using the generator `mg`. The momenta `p` are in units of GeV/c. The `cosθ` is the cosine of the zenith
angle (thus cosθ=1 means vertical, and 0 means horizontal).
"""
function generate(mg::CRMuonGenerator, N::Integer)
    boxid = [findfirst(x->x>r, mg.boxCDF) for r in rand(N)]

    cosθ = Array{Float64}(undef, N)
    logp = Array{Float64}(undef, N)
    dθ = mg.cosθlim[2]-mg.cosθlim[1]
    dlogp = mg.logPlim[2]-mg.logPlim[1]
    for i=1:N
        k = boxid[i]
        boxi = 1 + ((k-1) % mg.Np)
        boxj = 1 + div(k-1, mg.Np)
        B = 0.5*(mg.pij[k,1] + mg.pij[k,2])
        A = 1 - B
        Fy = rand()
        y = (sqrt(B^2 + 4A*Fy)-B)/2A
        cosθ[i] = y*dθ + mg.cosθlim[boxj]

        Bx = (1-y)*mg.pij[k,1] + y*mg.pij[k,3]
        Ax = 0.5*((1-y)*(mg.pij[k,2]-mg.pij[k,1]) + y*(mg.pij[k,4]-mg.pij[k,3]))
        rescale = Ax+Bx  # replace conditional prob with given y.
        Ax /= rescale
        Bx /= rescale
        Fx = rand()
        x = (sqrt(Bx^2 + 4Ax*Fx)-Bx)/2Ax
        logp[i] = x*dlogp + mg.logPlim[boxi]
    end
    p = exp.(logp)
    p, cosθ
end


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

        # Random position in the tube plane
        r = sqrt(rand())*tube_radius
        ψ = 2π*rand()
        x0 = r.*cos.(ψ)
        y0 = r.*sin.(ψ)

        # Where does the path cross the plane ⟂ the tube passing through the origin?
        origin = [cϕ*cθ*x0-sϕ*y0, sϕ*cθ*x0+cϕ*y0, -sθ*x0]

        line = Line(n, origin)
        L[i] = path(obj, line)
    end
    L
end
