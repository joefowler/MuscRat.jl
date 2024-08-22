using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
    MuscRat.Gamma,
)

function plot_all_spectra(alt=1462u"m")
    gp1 = MuscRat.OLDParmaGenerator(MuscRat.µplus)
    gm1 = MuscRat.OLDParmaGenerator(MuscRat.µminus)
    gp = ParmaGenerator(MuscRat.µplus; alt=alt)
    gm = ParmaGenerator(MuscRat.µminus; alt=alt)
    oldparma = sum(gp1.relativeprob, dims=2)*(gm1.flux*1u"s*cm^2" |> NoUnits) 
    oldparma += sum(gm1.relativeprob, dims=2)*(gm1.flux*1u"s*cm^2" |> NoUnits)
    parma = sum(gp.relativeprob, dims=2)*(gm.flux*1u"s*cm^2" |> NoUnits) 
    parma += sum(gm.relativeprob, dims=2)*(gm.flux*1u"s*cm^2" |> NoUnits)
    depth = 1u"g/cm^2"*MuscRat.Parma_depth(alt/1u"km" |> NoUnits)
    gr = CRMuonGenerator(100, 100; y=depth, Pmin=0.045u"GeV/c", useReyna=true)
    gc = CRMuonGenerator(100, 100; y=depth, Pmin=0.045u"GeV/c", useReyna=false)
    
    clf()
    P_GeV = exp.(0.5*(gp.logPGeVlim[2:end].+gp.logPGeVlim[1:end-1]))
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    # E_MeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"MeV" .|> NoUnits
    loglog(E_GeV, parma, "-", label="PARMA")
    loglog(E_GeV, oldparma, "-", label="PARMA OLD")

    P_GeV = exp.(gr.logPGeVlim[2:end])
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    reyna = sum(gr.relativeprob, dims=2)*(gr.flux*1u"s*cm^2" |> NoUnits)
    loglog(E_GeV, reyna, "-", label="Reyna")
    
    P_GeV = exp.(gc.logPGeVlim[2:end])
    E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
    # dE_MeV = 1e3*vcat(diff(E_GeV)..., E_GeV[end]-E_GeV[end-1])
    chatz = sum(gc.relativeprob, dims=2)*(gc.flux*1u"s*cm^2" |> NoUnits)
    loglog(E_GeV, chatz, "-", label="Chatzidakis")

    xlim([1e-4, 1e3])
    ylim([1e-9, 1e-3])
    grid()
    legend()
    xlabel("Energy (GeV)")
    ylabel("Spectrum (counts per cm\$^2\$ per s per LOGARITHMIC E BIN)")
    title("CR µ± spectra at $alt ASL or $depth")
end


function plot_muons_elevation()
    clf()
    for y0 in LinRange(880, 1030, 5)
        gc = CRMuonGenerator(100, 100; y=y0*u"g/cm^2", Pmin=0.045u"GeV/c", useReyna=false)
        P_GeV = exp.(gc.logPGeVlim[2:end])
        E_GeV = MuscRat.muon_E.(P_GeV*1u"GeV/c")/1u"GeV" .|> NoUnits
        chatz = sum(gc.relativeprob, dims=2)*(gc.flux*1u"s*cm^2" |> NoUnits)
        loglog(E_GeV, chatz, "-", label="Chatzidakis")
        @show sum(chatz)
    end
    xlabel("Energy (GeV)")
    ylabel("Flux (counts per s per cm^2)")
    title("Charzidakis µ± spectrum at depth [880...1030] g/cm^2")
    nothing
end

function plot_muons_elevation_parma(;nucl=false, electron=false, gamma=false)
    clf()
    lat, long, alt = 40, -105, 0.0
    E_GeV = 10 .^ LinRange(-2, 3, 101)
    particle1, particle2 = MuscRat.µminus, MuscRat.µplus
    if electron
        particle1, particle2 = MuscRat.Electron, MuscRat.Proton
        E_GeV = 10 .^ LinRange(-5, 3, 161)
    elseif nucl
        particle1, particle2 = MuscRat.Proton, MuscRat.Neutron
        E_GeV = 10 .^ LinRange(-6, 3, 181)
    elseif gamma
        particle1, particle2 = MuscRat.Gamma, nothing
        E_GeV = 10 .^ LinRange(-5, 3, 181)
    end
    E_MeV = E_GeV * 1e3
    dE_MeV = vcat(diff(E_MeV)..., E_MeV[end]-E_MeV[end-1])
    
    flux = Float64[]
    for alt in LinRange(1.46, 0, 5)
        obs = MuscRat.CRObserver(2022, 5, 23, lat, long, alt)
        s1 = MuscRat.CRspectrum(E_MeV, obs, particle1)
        if particle2 === nothing
            s2 = zero(s1)
        else
            s2 = MuscRat.CRspectrum(E_MeV, obs, particle2)
        end
        # Mult by dE_MeV to get per logarithmic bin, for comparison 
        spectrum = (s1 .+ s2) .* dE_MeV
        loglog(E_GeV, spectrum, "-")
        push!(flux, sum(spectrum))
    end
    xlabel("Energy (GeV)")
    ylabel("Flux (counts per s per cm^2 per log-bin)")
    flux
end
