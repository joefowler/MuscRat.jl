using PyPlot, Unitful, MuscRat, Statistics

particles = (
    MuscRat.µplus,
    MuscRat.µminus,
)

TKID = Box([5, 5, 1.5]*u"mm")
area = MuscRat.tube_area(TKID)
Ttotal = 5e7u"s"
Pmin = 10.0u"MeV/c"
loss_in_concrete = 0.0u"MeV"

pparma = Dict{MuscRat.Particle,Vector}()
cosθparma = Dict{MuscRat.Particle,Vector}()
for particle in particles
    g = ParmaGenerator(particle; Pmin=Pmin)
    rate = g.flux*area
    N = Int(round(Ttotal*rate))
    println("Generating $N $particle cosmic rays...")
    @time pparma[particle], cosθparma[particle] = generate(g, N)
end

p = vcat(pparma[MuscRat.µplus], pparma[MuscRat.µminus])
e = MuscRat.muon_E.(p) .- loss_in_concrete
e[e .< 0u"GeV"] .= 0u"GeV"
p = MuscRat.muon_p.(e)
allp = [p]
allcosθ = [vcat(cosθparma[MuscRat.µplus], cosθparma[MuscRat.µminus])]
algorithms = ["PARMA"]

Pmin = uconvert(u"GeV/c", Pmin)
generator = CRMuonGenerator(100, 100; Pmin=Pmin, useReyna=true);
rate = generator.flux*area
N = Int(round(Ttotal*rate))
println("Generating $N cosmic rays from Reyna...")
@time p,cosθ = generate(generator, N);
e = MuscRat.muon_E.(p) .- loss_in_concrete
e[e .< 0u"GeV"] .= 0u"GeV"
p = MuscRat.muon_p.(e)
push!(allp, p)
push!(allcosθ, cosθ)
push!(algorithms, "Reyna")

generator = CRMuonGenerator(100, 100; Pmin=Pmin, useReyna=false, y=1000.0u"g/cm^2");
rate = generator.flux*area
N = Int(round(Ttotal*rate))
println("Generating $N cosmic rays from Chatzidakis...")
@time p,cosθ = generate(generator, N);
e = MuscRat.muon_E.(p) .- loss_in_concrete
e[e .< 0u"GeV"] .= 0u"GeV"
p = MuscRat.muon_p.(e)
push!(allp, p)
push!(allcosθ, cosθ)
push!(algorithms, "Chatzidakis")

println("Computing path lengths")
@time total_paths = [MuscRat.path_values(TKID, cosθ) for cosθ in allcosθ]
println("Computing energy loss given paths")

eloss_tkid = Δprob(:Silicon, :µ)
@time Loss_tkid = [eloss_tkid.(allp[i], total_paths[i]) for i=1:3]

function plot_distributions(total_paths, Loss_tkid, lossMax=5.0)
    clf()
    ax1 = subplot(311); loglog(); xlabel("µ± momentum (GeV/c)")
    ylabel("Counts / GeV/c / cm\$^2\$ / second")
    ax2 = subplot(312); loglog(); xlabel("µ± kinetic energy (GeV)")
    ylabel("Counts / GeV / cm\$^2\$ / second")
    ax3 = subplot(313); semilogy(); xlabel("Energy loss in device (MeV)")
    ylabel("Counts / MeV / second")
    Nbins = 2000
    Pmax = 200.0
    Emax = 200.0
    lossbins = 500
    area = 0.5^2
    for i=1:3
        plen = total_paths[i]
        pGeV = allp[i]/1u"GeV/c"
        EGeV = MuscRat.muon_E.(allp[i])/1u"GeV"
        lossMeV = Loss_tkid[i][plen.>0u"cm"]/Unitful.MeV .|> NoUnits

        sca(ax1)
        binwidth = Pmax/Nbins
        w = ones(length(pGeV))/(Ttotal/1u"s")/binwidth/area
        hist(pGeV, Nbins, [0, Pmax], histtype="step", weights=w)

        sca(ax2)
        binwidth = Emax/Nbins
        w = ones(length(EGeV))/(Ttotal/1u"s")/binwidth/area
        hist(EGeV, Nbins, [0, Emax], histtype="step", weights=w)

        sca(ax3)
        binwidth = lossMax/lossbins
        w = ones(length(lossMeV))/(Ttotal/1u"s")/binwidth
        hist(lossMeV, lossbins, [0, lossMax], histtype="step", weights=w, label=algorithms[i])
    end
    legend()
end
plot_distributions(total_paths, Loss_tkid)