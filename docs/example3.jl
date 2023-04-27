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

loss = Dict{Symbol, Vector}()
meanloss = Dict{Symbol, Vector}()
for k in keys(solids)
    plen = total_paths[k]
    use = plen .> 0cm
    thisloss1 = probloss_tkid.(p[use], plen[use])
    loss[k] = uconvert.(MeV, thisloss1)

    thisloss2 = lossrate_tkid.(p[use]).*(total_paths[k])[use]
    @show thisloss1[1:5]
    @show thisloss2[1:5]
    meanloss[k] = uconvert.(MeV, thisloss2)
end

names = Dict(
    :thick_tkid => "1.5 mm TKID",
    :thin_tkid => "0.5 mm TKID",
)

objcolors = Dict(
    :thick_tkid => "red",
    :thin_tkid => "orange",
)

function plot_results(solids, loss, smear=0.0)
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
        lw = (k == :cylH2 ? 2 : 1)
        tp_cm = convert.(Float64, total_paths[k]/1cm)
        @show tp_cm[1:10]
        c, _, _ = hist(tp_cm[tp_cm.>0], 500, [0,convert(Float64,Pmax/1cm)], histtype="step", color=objcolors[k], label=names[k], lw=lw)

        sca(ax2)
        weight = 1/(total_time*Δbin)
        Ndetector_hits = length(loss[k])
        L = loss[k]
        if smear > 0
            L = L .* exp.(smear*randn(Ndetector_hits))
        end
        wtunits = 1u"1/s/MeV"
        c, _, _ = hist(L/MeV, N_lbins, [0,Lmax/MeV], histtype="step", weights=zeros(length(L)).+weight/wtunits,
                        color=objcolors[k], label=names[k], lw=lw)
        hitrate = Ndetector_hits/total_time
        @show k, hitrate

        Lm = meanloss[k]
        if smear > 0
            Lm = Lm .* exp.(smear*randn(Ndetector_hits))
        end
        c, _, _ = hist(Lm/MeV, N_lbins, [0,Lmax/MeV], histtype="step", weights=zeros(length(L)).+weight/wtunits,
                        color=objcolors[k], lw=lw/2, alpha=0.5)

        xlim([0, Lmax/Unitful.MeV])
    end
    legend()
    tight_layout()
    nothing
end

plot_results(solids, loss)

# using HDF5
# function store_results(solids, loss, smear=0.0; NaI=true)
#     if NaI
#         fname = "loss_spectra_NaI_10k.hdf5"
#         detectors = (:sphere, :cylH3, :cylV3, :cylH1, :cylH2, :cylV1, :cylV2)
#         detectors = (:cylH2, )
#         Lmax = 100 # MeV
#         N_lbins = 10000
#     else
#         fname = "loss_spectra_TKID.hdf5"
#         detectors = (:thick_tkid, :thin_tkid)
#         Lmax = 5 # MeV
#         N_lbins = 5000
#     end
#     figure(2)

#     h5open(fname, "w") do file
#         write(file, "Ebin_range", [0,Lmax])
#         for k in detectors
#             flux_in_tube = generator.flux*MuscRat.tube_area(solids[k]) # Units are µ per second
#             total_time = N/flux_in_tube
#             Δbin = Lmax/N_lbins # MeV

#             weight = 1/(total_time*Δbin)
#             Ndetector_hits = length(loss[k])
#             L = loss[k]
#             if smear > 0
#                 L = L .* exp.(smear*randn(Ndetector_hits))
#             end
#             c, _, _ = hist(L, N_lbins, [0,Lmax], histtype="step", weights=weight.+zero(L))
#             @show Ndetector_hits*weight*Δbin, sum(c)*Δbin
#             write(file, string(k), c)
#         end
#     end
#     close(2)
#     nothing
# end

