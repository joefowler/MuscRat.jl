using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
    MuscRat.Gamma,
)
weights = Dict(
    MuscRat.µplus => 1,
    MuscRat.µminus => 1,
    MuscRat.Electron => 1,
    MuscRat.Positron => 1,
    MuscRat.Gamma => 5,
)
Pmin = Dict(
    MuscRat.µplus => 1.0u"MeV/c",
    MuscRat.µminus => 1.0u"MeV/c",
    MuscRat.Electron => 1.0u"MeV/c",
    MuscRat.Positron => 1.0u"MeV/c",
    MuscRat.Gamma => 1.0u"MeV/c",
)

TKID = Box([5, 5, 1.5]*u"mm")
area = MuscRat.tube_area(TKID)
Ttotal = 1e7u"s"

p = Dict{MuscRat.Particle,Vector}()
cosθ = Dict{MuscRat.Particle,Vector}()
for particle in particles
    g = ParmaGenerator(particle; Pmin=Pmin[particle])
    rate = g.flux*area
    weight = weights[particle]
    N = Int(round(Ttotal*rate/weight))
    println("Generating $N $particle cosmic rays with weight $weight...")
    p[particle], cosθ[particle] = generate(g, N)
end

clf()
ax1 = subplot(211); loglog(); grid()
xlabel("Energy (GeV)")
ylabel("CR flux (counts per cm^2 per second per MeV)")
ax2 = subplot(212); semilogy(); grid()
xlabel("Energy")
ylabel("CR flux (counts per cm^2 per second per MeV)")
plotinfo = Dict(:emax=>500, :nbins=>20000, :nlogbins=>150, :logErange=>[-3, 3])

for (i, particle) in enumerate(particles)
    mass = MuscRat.masses[particle]
    EGeV = p[particle]/1u"GeV/c"
    if mass > 0u"MeV/c^2"
        EGeV = mass/1u"GeV/c^2"*(sqrt.(1 .+ (p[particle]/mass/Unitful.c) .^ 2) .- 1)
    end
    EGeV = uconvert.(NoUnits, EGeV)
    binwidth = 1000*plotinfo[:emax]/plotinfo[:nbins]
    w = ones(length(EGeV))*(1u"s"/Ttotal)*weights[particle]/binwidth/0.25
    sca(ax1)
    hist(EGeV, plotinfo[:nbins], [0, plotinfo[:emax]], histtype="step", weights=w, label=string(particle), color="C$i")

    # Read and display the PARMA model
    results = MuscRat.readParma(particle)
    model_E, fluxes = results[9:10]
    x = model_E/Unitful.GeV .|> NoUnits
    y = fluxes/1u"1/cm^2/s/MeV" .|> NoUnits
    plot(x, y, "--", color="C$i")
    
    sca(ax2)
    dE = diff(plotinfo[:logErange])[1]
    bfrac = dE/plotinfo[:nlogbins]
    EMeV = 1000EGeV
    logbinwidth = EMeV*log(10.0)*bfrac
    w .*= binwidth./logbinwidth

    hist(log10.(EGeV), plotinfo[:nlogbins], plotinfo[:logErange], histtype="step", weights=w, label=string(particle), color="C$i")
    # Ddisplay the PARMA model
    plot(log10.(x), y, "--", color="C$i")
    
end
sca(ax1)
xlim([.01, 500])
legend()
sca(ax2)
x = LinRange(-3, 3, 7)
xt = ["1 MeV", "10 MeV", "100 MeV", "1 GeV", "10 GeV", "100 GeV", "1 TeV"]
xticks(x, xt)
xlim([-3, 3])