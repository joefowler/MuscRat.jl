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
    MuscRat.Electron => 30,
    MuscRat.Positron => 10,
    MuscRat.Gamma => 3000,
)
Pmin = Dict(
    MuscRat.µplus => 10.0u"MeV/c",
    MuscRat.µminus => 10.0u"MeV/c",
    MuscRat.Electron => 100.0u"MeV/c",
    MuscRat.Positron => 100.0u"MeV/c",
    MuscRat.Gamma => 300.0u"MeV/c",
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

clf(); loglog()
xlabel("Energy (GeV)")
ylabel("CR flux (counts per cm^2 per second per MeV)")
for particle in particles
    mass = MuscRat.masses[particle]
    EGeV = p[particle]/1u"GeV/c"
    if mass > 0u"MeV/c^2"
        EGeV = mass/1u"GeV/c^2"*(sqrt.(1 .+ (p[particle]/mass/Unitful.c) .^ 2) .- 1)
    end
    EGeV = uconvert.(NoUnits, EGeV)
    binwidth = 100 # MeV
    w = ones(length(EGeV))*(1u"s"/Ttotal)*weights[particle]/binwidth/0.25
    hist(EGeV, 1000, [0, 100], histtype="step", weights=w, label=string(particle))
end
legend()
