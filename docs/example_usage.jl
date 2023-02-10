using PyPlot, MuscRat, Statistics

# Make 2 boxes: 5x5x1.5 mm and 5x5x0.5 mm
b3 = Box(.5, .5, .15)
b1 = Box(.5, .5, .05)

# Make 3 "NaI" objects: a cylinder of the true shape, and 2 boxes of same volume
NaI = HCylinder(3.81, 7.62)
v = volume(NaI)
w = sqrt(v/7.62)
NaBtall = Box(w, w, 7.62)
NaBfat = Box(7.62, 7.62, v/7.62^2)

solids = Dict(
    :thick_tkid => b3,
    :thin_tkid => b1,
    :cylinder => NaI,
    :tall_scint => NaBtall,
    :fat_scint => NaBfat
)

N = 1000000
generator = CRMuonGenerator(100, 100);
pGeV,cosθ = generate(generator, N);
pMeV = 1000pGeV

total_paths = Dict([(k,MuscRat.path_values(obj, cosθ)) for (k,obj) in solids])

eloss_tkid = Eloss_function(:Silicon, :µ)
eloss_scint = Eloss_function(:NaI, :µ)

LossRate_tkid = eloss_tkid.(pMeV)
LossRate_scint = eloss_scint.(pMeV)
lossrate = Dict([(k,LossRate_scint) for k in keys(solids)])
lossrate[:thin_tkid] = LossRate_tkid
lossrate[:thick_tkid] = LossRate_tkid

loss = Dict{Symbol, Vector{Float64}}()
for k in keys(solids)
    P = total_paths[k]
    loss[k] = (P.*lossrate[k])[P.>0]
end

names = Dict(
    :thick_tkid => "1.5 mm TKID",
    :thin_tkid => "0.5 mm TKID",
    :cylinder => "NaI cylinder",
    :tall_scint => "NaI tall box",
    :fat_scint => "NaI fat box"
)

function plot_NaI_results(solids, loss)
    nbins = 500
    brange = [0,100]
    db = (brange[2]-brange[1])/nbins
    clf()
    for k in (:cylinder, :tall_scint, :fat_scint)
        weight = db*generator.flux*MuscRat.tube_area(solids[k])/N
        c, _, _ = hist(loss[k], 500, [0,100], histtype="step", weights=weight.+zero(loss[k]), label=names[k])
        @show length(loss[k])*weight/db, sum(c)*db
    end
    title("Bethe-Bloch µ± loss in NaI of 3 volumes, equal shape")
    xlabel("Energy lost in scintillator (MeV)")
    ylabel("Counts per MeV per second")
    legend()
    tight_layout()
    nothing
end

plot_NaI_results(solids, loss)
