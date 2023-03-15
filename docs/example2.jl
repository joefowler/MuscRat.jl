using PyPlot, MuscRat, Statistics

particles = (
    MuscRat.Gamma,
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
)

TKID = Box([.5, .5, .05])
area = MuscRat.tube_area(TKID)*1u"cm^2"
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
