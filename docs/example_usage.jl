using PyPlot, MuscRat, Statistics, Unitful
import Unitful.cm, Unitful.mm, Unitful.MeV

# Make 7 "NaI" objects of equal volume: a sphere, three cylinders oriented vertically then horizontally
r = 3.81cm
h = 2r
NaIH2 = HCylinder(r, h)
NaIV2 = VCylinder(r, h)
NaIH1 = HCylinder(r/sqrt(2), 2h)
NaIV1 = VCylinder(r/sqrt(2), 2h)
NaIH3 = HCylinder(sqrt(2)*r, h/2)
NaIV3 = VCylinder(sqrt(2)*r, h/2)
v = volume(NaIH2)
rsphere = (3/(4π)*v)^(1/3)
NaISphere = Sphere(rsphere)

solids = Dict(
    :thick_tkid => Box(5mm, 5mm, 1.5mm),
    :thin_tkid => Box(5mm, 5mm, 0.5mm),
    :cylH1 => NaIH1,
    :cylH2 => NaIH2,
    :cylH3 => NaIH3,
    :cylV1 => NaIV1,
    :cylV2 => NaIV2,
    :cylV3 => NaIV3,
    :sphere => NaISphere,
)

N = 1000000
generator = CRMuonGenerator(100, 100);
@time p,cosθ = generate(generator, N);
println("Generated $N muons")

total_paths = Dict([(k,MuscRat.path_values(obj, cosθ)) for (k,obj) in solids])

eloss_tkid,_ = Eloss_functions(:Silicon, :µ)
eloss_scint,_ = Eloss_functions(:NaI, :µ)

LossRate_tkid = eloss_tkid.(p)
LossRate_scint = eloss_scint.(p)
lossrate = Dict([(k,LossRate_scint) for k in keys(solids)])
lossrate[:thin_tkid] = LossRate_tkid
lossrate[:thick_tkid] = LossRate_tkid

loss = Dict{Symbol, Vector}()
for k in keys(solids)
    plen = total_paths[k]
    loss[k] = uconvert.(MeV, plen.*lossrate[k])[plen.>0cm]
end

names = Dict(
    :thick_tkid => "1.5 mm TKID",
    :thin_tkid => "0.5 mm TKID",
    :sphere => "NaI sphere",
    :cylH1 => "NaI long horizontal cylinder",
    :cylH2 => "NaI equal horizontal cylinder",
    :cylH3 => "NaI short horizontal cylinder",
    :cylV1 => "NaI tall vertical cylinder",
    :cylV2 => "NaI equal vertical cylinder",
    :cylV3 => "NaI short vertical cylinder",
)

objcolors = Dict(
    :thick_tkid => "red",
    :thin_tkid => "orange",
    :sphere => "brown",
    :cylH1 => "navy",
    :cylH2 => "black",
    :cylH3 => "blue",
    :cylV1 => "darkturquoise",
    :cylV2 => "green",
    :cylV3 => "limegreen",
)

function plot_results(solids, loss, smear=0.0; NaI=true)
    clf()
    ax1 = subplot(211)
    if NaI
        title("Geometric paths for CR µ± through NaI of 4 shapes of equal volume")
        xlabel("Path through scintillator (cm)")
    else
        title("Geometric paths for CR µ± through TKID of 2 thicknesses")
        xlabel("Path through silicon TKID (cm)")
    end
    ax2 = subplot(212)
    if NaI
        title("Bethe-Bloch CR µ± loss in NaI of 4 shapes of equal volume")
        xlabel("Energy lost in scintillator (MeV)")
    else
        title("Bethe-Bloch CR µ± loss in TKID of 2 thicknesses")
        xlabel("Energy lost in silicon TKID (MeV)")
    end
    ylabel("Counts per MeV per second")

    if NaI
        detectors = (:cylH1, :cylH2, :cylH3, :cylV1, :cylV2, :cylV3, :sphere)
        Pmax, Lmax = 2smallest_radius(solids[:cylH1]), 150MeV
    else
        detectors = (:thick_tkid, :thin_tkid)
        Pmax, Lmax = 7.5mm, 5MeV
    end

    for k in detectors
        flux_in_tube = generator.flux*MuscRat.tube_area(solids[k]) # Units are µ per second
        total_time = N/flux_in_tube
        N_lbins = 500
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
    end
    legend()
    tight_layout()
    nothing
end

plot_results(solids, loss; NaI=false)

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

