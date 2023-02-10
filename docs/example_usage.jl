using PyPlot, MuscRat, Statistics

# Make 2 boxes: 5x5x1.5 mm and 5x5x0.5 mm
b3 = Box(.5, .5, .15)
b1 = Box(.5, .5, .05)

# Make 4 "NaI" objects: a cylinder of the true shape, plus a sphere and 2 boxes of same volume
NaI = HCylinder(3.81, 7.62)
v = volume(NaI)
NaISphere = Sphere((3/(4π)*v)^(1/3))
w = sqrt(v/7.62)
NaBtall = Box(w, w, 7.62)
NaBfat = Box(7.62, 7.62, v/7.62^2)

solids = Dict(
    :thick_tkid => b3,
    :thin_tkid => b1,
    :cylinder => NaI,
    :sphere => NaISphere,
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
    :sphere => "NaI sphere",
    :tall_scint => "NaI tall box",
    :fat_scint => "NaI fat box"
)

objcolors = Dict(
    :thick_tkid => "red",
    :thin_tkid => "orange",
    :cylinder => "black",
    :sphere => "brown",
    :tall_scint => "blue",
    :fat_scint => "cyan"
)

function plot_NaI_results(solids, loss, smear=0.0)
    nbins = 500
    brange = [0,100]
    db = (brange[2]-brange[1])/nbins
    clf()
    ax1 = subplot(211)
    title("Geometric paths for CR µ± through NaI of 4 shapes of equal volume")
    xlabel("Path through scintillator (cm)")
    ax2 = subplot(212)
    title("Bethe-Bloch CR µ± loss in NaI of 4 shapes of equal volume")
    xlabel("Energy lost in scintillator (MeV)")
    ylabel("Counts per MeV per second")

    for k in (:cylinder, :sphere, :tall_scint, :fat_scint)
        sca(ax1)
        tp = total_paths[k]
        lw = k == :cylinder ? 2 : 1
        c, _, _ = hist(tp[tp.>0], 600, [0,12], histtype="step", color=objcolors[k], label=names[k], lw=lw)
        sca(ax2)
        weight = generator.flux*MuscRat.tube_area(solids[k])/(N*db)
        Nv = length(loss[k])
        L = loss[k] .* exp.(smear*randn(Nv))
        c, _, _ = hist(L, 500, [0,100], histtype="step", weights=weight.+zero(loss[k]), 
                        color=objcolors[k], label=names[k], lw=lw)
        @show Nv*weight*db, sum(c)*db
    end
    legend()
    tight_layout()
    nothing
end

plot_NaI_results(solids, loss)
