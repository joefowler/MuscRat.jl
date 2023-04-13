using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
    MuscRat.Gamma,
)

TKID = Box([5, 5, 1.5]*u"mm")
area = MuscRat.tube_area(TKID)
Ttotal = 1e8u"s"
Pmin = 10.0u"MeV/c"

p = Dict{MuscRat.Particle,Vector}()
cosθ = Dict{MuscRat.Particle,Vector}()
for particle in particles
    g = ParmaGenerator(particle; Pmin=Pmin)
    rate = g.flux*area
    N = Int(round(Ttotal*rate))
    println("Generating $N $particle cosmic rays...")
    p[particle], cosθ[particle] = generate(g, N)
end
clf(); loglog()
xlabel("Energy (GeV)")
for particle in particles
    mass = MuscRat.masses[particle]
    EGeV = p[particle]/1u"GeV/c"
    if mass > 0u"MeV/c^2"
        EGeV = mass/1u"GeV/c^2"*(sqrt.(1 .+ (p[particle]/mass/Unitful.c) .^ 2) .- 1)
    end
    EGeV = uconvert.(NoUnits, EGeV)
    hist(EGeV, 1000, [0, 100], histtype="step")
end
