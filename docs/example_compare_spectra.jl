using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
    MuscRat.Gamma,
)

function plot_all_spectra()
    gp = ParmaGenerator(MuscRat.µplus)
    gm = ParmaGenerator(MuscRat.µminus)
    # gp = ParmaGenerator(MuscRat.µplus; Pmin=0.045u"GeV/c")
    # gm = ParmaGenerator(MuscRat.µminus; Pmin=0.045u"GeV/c")
    parma = sum(gp.relativeprob, dims=2)*(gm.flux*1u"s*cm^2" |> NoUnits) 
    parma += sum(gm.relativeprob, dims=2)*(gm.flux*1u"s*cm^2" |> NoUnits)
    gr = CRMuonGenerator(100, 100; y=1030u"g/cm^2", Pmin=0.045u"GeV/c", useReyna=true)
    gc = CRMuonGenerator(100, 100; y=1030u"g/cm^2", Pmin=0.045u"GeV/c", useReyna=false)
    reyna = sum(gr.relativeprob, dims=2)*(gr.flux*1u"s*cm^2" |> NoUnits)
    chatz = sum(gc.relativeprob, dims=2)*(gc.flux*1u"s*cm^2" |> NoUnits)
    
    clf()
    P_GeV = exp.(0.5*(gp.logPGeVlim[2:end].+gp.logPGeVlim[1:end-1]))
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    # E_MeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"MeV" .|> NoUnits
    loglog(E_GeV, parma, "-", label="PARMA")

    P_GeV = exp.(gr.logPGeVlim[2:end])
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    loglog(E_GeV, reyna, "-", label="Reyna")

    P_GeV = exp.(gc.logPGeVlim[2:end])
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    loglog(E_GeV, chatz, "-", label="Chatzidakis")

    xlim([1e-4, 1e3])
    ylim([1e-9, 1e-3])
    grid()
    legend()
    xlabel("Energy (GeV)")
    ylabel("Spectrum (counts per cm\$^2\$ per s per LOGARITHMIC E BIN)")
end


function plot_muons_elevation()
    clf()

    for y0 in LinRange(850, 1030, 5)
        gc = CRMuonGenerator(100, 100; y=y0*u"g/cm^2", Pmin=0.045u"GeV/c", useReyna=false)
        P_GeV = exp.(gc.logPGeVlim[2:end])
        E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
        chatz = sum(gc.relativeprob, dims=2)*(gc.flux*1u"s*cm^2" |> NoUnits)
        loglog(E_GeV, chatz, "-", label="Chatzidakis")
        @show sum(chatz)
    end
    xlabel("Energy (GeV)")
    ylabel("Flux (counts per s per cm^2)")
end
