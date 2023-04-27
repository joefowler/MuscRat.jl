using PyPlot, MuscRat, Statistics, Unitful
import Unitful.cm, Unitful.mm, Unitful.MeV

solids = Dict(
    :thick_tkid => Box(5mm, 5mm, 1.5mm),
    :thin_tkid => Box(5mm, 5mm, 0.5mm),
)

N = 1000000
generator = CRMuonGenerator(100, 100);
@time p,cosθ = generate(generator, N);
println("Generated $N muons")

total_paths = Dict([(k,MuscRat.path_values(obj, cosθ)) for (k,obj) in solids])

lossrate_tkid, probloss_tkid = Eloss_functions(:Silicon, :µ)
loss_distribution = MuscRat.Eloss_µSi()

loss = Dict{Symbol, Vector}()
mploss = Dict{Symbol, Vector}()
meanloss = Dict{Symbol, Vector}()
for k in keys(solids)
    plen = total_paths[k]
    use = plen .> 0cm
    thisloss1 = probloss_tkid.(p[use], plen[use])
    mploss[k] = uconvert.(MeV, thisloss1)

    thisloss2 = lossrate_tkid.(p[use]).*plen[use]
    meanloss[k] = uconvert.(MeV, thisloss2)

    thisloss3 = loss_distribution.(p[use], plen[use])
    loss[k] = uconvert.(MeV, thisloss3)
end

names = Dict(
    :thick_tkid => "1.5 mm TKID",
    :thin_tkid => "0.5 mm TKID",
)

objcolors = Dict(
    :thick_tkid => "red",
    :thin_tkid => "orange",
)

function plot_results(solids, loss, mploss, meanloss, smear=0.0)
    clf()
    ax1 = subplot(211)
    title("Geometric paths for CR µ± through TKID of 2 thicknesses")
    xlabel("Path through silicon TKID (cm)")
    ax2 = subplot(212)
    title("Bethe-Bloch CR µ± loss in TKID of 2 thicknesses")
    xlabel("Energy lost in silicon TKID (MeV)")
    ylabel("Counts per MeV per second")

    detectors = (:thick_tkid, :thin_tkid)
    Pmax, Lmax = 7.5mm, 5MeV

    for k in detectors
        flux_in_tube = generator.flux*MuscRat.tube_area(solids[k]) # Units are µ per second
        total_time = N/flux_in_tube
        N_lbins = 1000
        Δbin = Lmax/N_lbins

        sca(ax1)
        lw = 1
        tp_cm = convert.(Float64, total_paths[k]/1cm)
        @show tp_cm[1:10]
        c, _, _ = hist(tp_cm[tp_cm.>0], 500, [0,convert(Float64,Pmax/1cm)], histtype="step", color=objcolors[k], label=names[k], lw=lw)

        sca(ax2)
        weight = 1/(total_time*Δbin)
        Ndetector_hits = length(loss[k])

        L = loss[k]/MeV .|> NoUnits
        wtunits = 1u"1/s/MeV"
        c, _, _ = hist(L, N_lbins, [0,Lmax/MeV], histtype="step", weights=zeros(length(L)).+weight/wtunits,
                        color=objcolors[k], label=names[k], lw=2)
        hitrate = Ndetector_hits/total_time
        @show k, hitrate

        Lm = mploss[k]/MeV .|> NoUnits
        if smear > 0
            Lm = Lm .* exp.(smear*randn(Ndetector_hits))
        end
        c, _, _ = hist(Lm, N_lbins, [0,Lmax/MeV], histtype="step", weights=zeros(length(L)).+weight/wtunits,
                        color=objcolors[k], lw=lw/2, alpha=0.85, label="Most prob, smeared")

        Lm = meanloss[k]/MeV .|> NoUnits
        if smear > 0
            Lm = Lm .* exp.(smear*randn(Ndetector_hits))
        end
        c, _, _ = hist(Lm, N_lbins, [0,Lmax/MeV], histtype="step", weights=zeros(length(L)).+weight/wtunits,
                        color=objcolors[k], lw=lw/2, alpha=0.5, label="Bethe-Bloch, smeared")

        xlim([0, Lmax/Unitful.MeV])
    end
    semilogy()
    legend()
    tight_layout()
    nothing
end

plot_results(solids, loss, mploss, meanloss, 0.15)

# using HDF5
# function store_results(solids, loss)
#     fname = "loss_spectra_TKID.hdf5"
#     detectors = (:thick_tkid, :thin_tkid)
#     Lmax = 5 # MeV
#     N_lbins = 5000
#     figure(2)

#     h5open(fname, "w") do file
#         write(file, "Ebin_range", [0,Lmax])
#         for k in detectors
#             flux_in_tube = generator.flux*MuscRat.tube_area(solids[k]) # Units are µ per second
#             flux_in_tube = uconvert(u"s^-1", flux_in_tube)
#             total_time = N/flux_in_tube
#             Δbin = Lmax/N_lbins # MeV
#             @show flux_in_tube, total_time, Δbin

#             weight = 1u"s"/(total_time*Δbin)
#             Ndetector_hits = length(loss[k])
#             L = loss[k]/1u"MeV" .|> NoUnits
#             L = L[L .> 0]
#             c, _, _ = hist(L, N_lbins, [0,Lmax], histtype="step", weights=weight.+zero(L))
#             @show Ndetector_hits*weight*Δbin, sum(c)*Δbin
#             write(file, string(k), c)
#         end
#     end
#     nothing
# end

# store_results(solids, loss)
