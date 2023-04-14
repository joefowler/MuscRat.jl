using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
    MuscRat.Gamma,
)

function plot_parma_spectra()
    clf()
    fluxes = Dict()
    for (i, p) in enumerate(particles)
        Np, Nang, Pmin, Pmax, logPGeVlim, mass, cosθlim, spectrum, KE, flux = MuscRat.readParma(p)
        fluxes[p] = flux
        KEMeV = KE/1u"MeV" .|> NoUnits
        pMeV = exp.(logPGeVlim)*1000
        fluxn = flux/1u"1/cm^2/s/MeV" .|> NoUnits
        loglog(KEMeV, fluxn, label=string(p), color="C$i")

        angle_avg = 2π*Unitful.sr*mean(spectrum, dims=2)[:,1]/1u"1/cm^2/s/MeV" .|> NoUnits
        loglog(KEMeV, angle_avg, "--", color="C$i")
        if p == MuscRat.µminus
            f = (fluxes[MuscRat.µplus]+fluxes[MuscRat.µminus])/1u"1/cm^2/s/MeV" .|> NoUnits
            loglog(KEMeV, f, color="k", label="All µ")
        end
    end
    xlim([1, 1e5])
    ylim([1e-9, 1e-2])
    grid()
    legend()
    xlabel("Energy (MeV)")
    ylabel("Spectrum (counts per cm\$^2\$ per s per MeV)")
end