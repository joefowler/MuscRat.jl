using PyPlot, MuscRat, Statistics

particles = (
    MuscRat.Gamma,
    MuscRat.µplus,
    MuscRat.µminus,
    MuscRat.Electron,
    MuscRat.Positron,
)

for particle in particles
    g = ParmaGenerator(particle)
    @show g.flux
end
