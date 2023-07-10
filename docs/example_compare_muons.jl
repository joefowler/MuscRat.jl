using PyPlot, Unitful, MuscRat, Statistics

Ttotal = 1e9u"s"

function run_simulation()
    TKID = Box([5, 5, 1.5]*u"mm")
    area = MuscRat.tube_area(TKID)
    Pmin = 10.0u"MeV/c"
    loss_in_concrete = 0.0u"MeV"

    particles = (
        MuscRat.µplus,
        MuscRat.µminus,
    )
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

    for depth in (1000u"g/cm^2", 850u"g/cm^2")
        generator = CRMuonGenerator(100, 100; Pmin=Pmin, useReyna=false, y=depth);
        rate = generator.flux*area
        N = Int(round(Ttotal*rate))
        println("Generating $N cosmic rays from Chatzidakis at depth $depth...")
        @time p,cosθ = generate(generator, N);
        e = MuscRat.muon_E.(p) .- loss_in_concrete
        e[e .< 0u"GeV"] .= 0u"GeV"
        p = MuscRat.muon_p.(e)
        push!(allp, p)
        push!(allcosθ, cosθ)
        push!(algorithms, "Chatzidakis $depth")
    end

    println("Computing path lengths")
    @time total_paths = [MuscRat.path_values(TKID, cosθ) for cosθ in allcosθ]
    println("Computing energy loss given paths")

    loss_distribution = MuscRat.Eloss_µSi()
    @time Loss_tkid = [loss_distribution.(allp[i], total_paths[i]) for i=1:4]
    allp, allcosθ, algorithms, total_paths, Loss_tkid
end
allp, allcosθ, algorithms, total_paths, Loss_tkid = run_simulation()

function plot_distributions(allp, total_paths, Loss_tkid, lossMax=5.0)
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
    for i=1:4
        plen = total_paths[i]
        pGeV = allp[i]/1u"GeV/c"
        EGeV = MuscRat.muon_E.(allp[i])/1u"GeV"
        lossMeV = Loss_tkid[i][plen.>0u"cm"]/Unitful.MeV .|> NoUnits
        # lossMeV .*= exp.(randn(length(lossMeV))*0.15)

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

plot_distributions(allp, total_paths, Loss_tkid)

function save_distributions(allp, allcosθ, total_paths, Loss_tkid, algorithms)
    h5open("muon_CR_loss_4models.hdf5", "w") do h
        for (i, k) in enumerate(algorithms)
            name = split(k, " g ")[1]
            name = replace(name, " "=>"_")
            g = create_group(h, name)

            use = total_paths[i] .> 0u"mm"
            write(g, "p_GeV", allp[i][use]/1u"GeV/c")
            write(g, "costheta", allcosθ[i][use])
            write(g, "path_mm", total_paths[i][use]/1u"mm")
            write(g, "loss_keV", Loss_tkid[i][use]/1u"keV")
        end
    end
end
function save_histograms(total_paths, Loss_tkid, algorithms)
    figure(2)
    clf()
    lossmax = 6000.0
    lossbins = 600
    binwidth = lossmax/lossbins/1000.0
    h5open("muon_CR_loss_histogram_4models.hdf5", "w") do h
        for (i, k) in enumerate(algorithms)
            name = split(k, " g ")[1]
            name = replace(name, " "=>"_")

            use = total_paths[i] .> 0u"mm"
            w = (Ttotal/1u"s")*binwidth
            c, b, _ = hist(Loss_tkid[i][use]/1u"keV", lossbins, [0,lossmax], histtype="step")
            write(h, name, c/w)
            ds = h[name]
            attr = Dict("binwidth"=>10.0, "Emax"=>lossmax, "Eunits"=>"keV", "histunits"=>"counts/sec/MeV")
            for (k,v) in attr
                attributes(ds)[k]=v
            end
        end
    end
end
save_histograms(total_paths, Loss_tkid, algorithms)
Ttotal <= 3e7u"s" && save_distributions(allp, allcosθ, total_paths, Loss_tkid, algorithms)
