using MuscRat, PyPlot, Quadrature, FastGaussQuadrature
using LinearAlgebra

function flux_vs_depth(obs::MuscRat.CRObserver, id::MuscRat.ParmaParticle)
    alt = LinRange(0, 10, 151)
    depth = MuscRat.Parma_depth.(alt)
    f(x,p) = MuscRat.CRspectrum(x, p[1],obs, id)
    flux = Float64[]
    Emin, Emax = 1e-2, 1e7
    if id == MuscRat.ppGamma
        Emin = 1e-4
    end
    for (i, d) in enumerate(depth)
        prob = QuadratureProblem(f, Emin, Emax, [d])
        result = solve(prob, HCubatureJL(), reltol=1e-3)
        push!(flux, result.u)
        if id == MuscRat.ppGamma
            flux[end] += MuscRat.Parma_get511flux(obs.Windex, obs.cutoffRigidity, d)
        end
    end
    alt, depth, flux
end

function dose_vs_depth(obs::MuscRat.CRObserver, id::MuscRat.ParmaParticle)
    alt = LinRange(0, 2, 31)
    depth = MuscRat.Parma_depth.(alt)
    f(x,p) = MuscRat.CRspectrum(x, p[1],obs, id)*x
    flux = Float64[]
    Emin, Emax = 1e0, 1e7
    if id == MuscRat.ppGamma
        Emin = 1e-4
    end
    for (i, d) in enumerate(depth)
        prob = QuadratureProblem(f, Emin, Emax, [d])
        result = solve(prob, HCubatureJL(), reltol=1e-3)
        push!(flux, result.u)
    end
    alt, depth, flux
end

function all_fluxes(obs::MuscRat.CRObserver, normalize::Bool=false, dose::Bool=false)
    fvdepth = flux_vs_depth
    if dose
        fvdepth = dose_vs_depth
    end
    a,d,fluxG = fvdepth(obs, MuscRat.ppGamma)
    a,d,fe = fvdepth(obs, MuscRat.ppElectron)
    a,d,fp = fvdepth(obs, MuscRat.ppPositron)
    a,d,fluxP = fvdepth(obs, MuscRat.ppProton)
    a,d,fluxN = fvdepth(obs, MuscRat.ppNeutron)
    a,d,fmp = fvdepth(obs, MuscRat.ppµplus)
    a,d,fmm = fvdepth(obs, MuscRat.ppµminus)
    for f in (fmm, fmp, fe, fp, fluxG, fluxP, fluxN)
        print(f[1])
    end
    fluxM = fmm .+ fmp
    if normalize
        fluxG ./= fluxG[1]
        fe ./= fe[1]
        fp ./= fp[1]
        fluxP ./= fluxP[1]
        fluxN ./= fluxN[1]
        fluxM ./= fluxM[1]
    end
    clf()
    plot(a, fluxM, ".-", label="µ±")
    plot(a, fe, ".-", label="e-")
    plot(a, fp, ".-", label="e+")
    plot(a, fluxG, ".-", label="γ")
    plot(a, fluxP, ".-", label="p")
    plot(a, fluxN, ".-", label="n")
    xlabel("Altitude (km a.s.l.)")
    if normalize
        if dose
            title("Cosmic-ray dose (angle-integrated) by particle, relative to sea level")
            ylabel("Total particle dose relative to sea level")
        else
            title("Cosmic-ray flux (angle-integrated) by particle, relative to sea level")
            ylabel("Total particle flux relative to sea level")
        end
    else
        if dose
            ylabel("Total particle dose (MeV cm\$^{-2}\$ s\$^{-1}\$)")
            title("Cosmic-ray dose (angle-integrated) by particle")
        else
            ylabel("Total particle flux (cm\$^{-2}\$ s\$^{-1}\$)")
            title("Cosmic-ray flux (angle-integrated) by particle")
        end
        semilogy()
    end
    legend()
end


function plot_spectrum(obs::MuscRat.CRObserver, depth::Real)
    E = 10 .^ LinRange(-2, 7, 181)
    clf()
    ids = (
        MuscRat.ppµminus,
        MuscRat.ppµplus,
        MuscRat.ppElectron,
        MuscRat.ppPositron,
        MuscRat.ppGamma,
        MuscRat.ppProton,
        MuscRat.ppNeutron,
        )
    labels = ("µ-", "µ+", "e-", "e+", "γ", "p", "n")
    for (id, label) in zip(ids, labels)
        s = [MuscRat.CRspectrum(x, obs, id, depth=depth) for x in E]
        if id == MuscRat.ppGamma
            bin = argmax(E .> 0.511) - 1
            binwid = E[bin+1]-E[bin]
            s[bin] += MuscRat.Parma_get511flux(obs.Windex, obs.cutoffRigidity, depth)/binwid
        end
        plot(E, s, label=label)
    end
    loglog()
    legend()
    xlabel("Energy (MeV)")
    ylabel("Flux (cm\$^{-2}\$⋅s\$^{-1}\$)⋅MeV\$^{-1}\$")
    ylim([1e-16, 1])
    xlim([1e-2, 1e7])
    title("Cosmic-ray spectrum at depth $(depth) g/cm\$^2\$")
end


function compare_angle(obs::MuscRat.CRObserver)
    cosθ, w = gausslegendre(65)
    id = MuscRat.ppProton
    for E in 10 .^ LinRange(0, 6, 7)
        s1 = MuscRat.CRspectrum(E, obs, id)
        s2points = MuscRat.CRangularSpectrum(E, cosθ, obs, id)
        s2 = 2π * dot(s2points, w)
        @show s1, s2/s1
    end
end