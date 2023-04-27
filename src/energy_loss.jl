import Unitful.MeV, Unitful.eV
using Interpolations

"""
    Eloss_functions(material=:Silicon, particle=:µ)

Returns two _functions_ `dEdx, Δprob`
1. `dEdX(p)` that accepts a momentum (MeV/c units) and returns the mean energy loss per unit length.
2. `Δprob(p, t)` that accepts a momentum and a material thickness `t` and returns the most probable loss.
The mean energy is according to the Bethe-Bloch formula (Eq 34.5 in the 2022 _Review of
Particle Properties_), while the most probable energy loss is Eq 34.12.

* `material` is the absorbing material. Allowed values: (:Silicon, :Copper, :NaI)
* `particle` is the particle being stopped. Allowed values: (:µ, :e, :p)
"""
function Eloss_functions(material=:Silicon, particle=:µ)
    massenergies = Dict(
        :µ => 105.65837MeV,
        :π => 139.5706MeV,
        :e => 0.5109989MeV,
        :p => 938.27208MeV,
    )
    Mc2 = massenergies[particle]
    K = 0.307075u"MeV/mol*cm^2"
    if material == :Silicon
        I = 175eV  # silicon excitation energy (eV, not MeV, for some reason!) See Fig 33.5
        ZoverA = 14/28.085u"g/mol"  # silicon charge-to-mass ratio
        ρ_absorber = 2.329u"g/cm^3"
    elseif material == :Copper
        I = 29*11eV
        ZoverA = 0.45636u"mol/g"
        ρ_absorber = 8.960u"g/cm^3"
    elseif material == :NaI
        I = .5*(140eV+500eV)  # roughly 140 eV and 500 eV for Na and I separately...so average them?
        ZoverA = 0.42697u"mol/g"
        ρ_absorber = 3.667u"g/cm^3"
    else
        error("material must be one of (:Silicon, :Copper, :NaI)")
    end

    # This is the Bethe-Bloch dE/dx value in units of MeV/cm
    dEdx_bethe = function(p::Unitful.Momentum)
        βγ = p*Unitful.c/Mc2
        γ = sqrt(1+βγ^2)
        β = βγ/γ
        mratio = massenergies[:e]/Mc2
        Wmax_factor = 1/sqrt(1+2γ*mratio+mratio^2)
        bracket = log(2massenergies[:e]*βγ^2*Wmax_factor/I)-β^2
        dedx = ρ_absorber*K*(ZoverA)*bracket/β^2
        uconvert(u"MeV/cm", dedx)
    end

    Δprob = function (p::Unitful.Momentum, thickness::Unitful.Length)
        if thickness ≤ 0u"mm"
            return 0
        end
        βγ = p*Unitful.c/Mc2
        γ = sqrt(1+βγ^2)
        β = βγ/γ
        ξ = 0.5*K*ZoverA*ρ_absorber*thickness/β^2
        j = 0.200
        bracket = log(2massenergies[:e]*βγ^2/I) + log(ξ/I) + j - β^2
        Δp = ξ*bracket
        uconvert(u"MeV", Δp)
    end
    return dEdx_bethe, Δprob
end

dEdx(material=:Silicon, particle=:µ)=Eloss_functions(material, particle)[1]
Δprob(material=:Silicon, particle=:µ)=Eloss_functions(material, particle)[2]

"""
Eloss_µSi()

Return a function f(e, path) to compute the random energy loss of a muon with energy `e` passing through a 
path of length `path` in silicon.

This function uses parameterized GEANT4 results taken at e=1 GeV and paths of 0.5 and 1.5 mm. We assume
the straggle function scales with Δprob(e, path), the most probable energy loss at energy `e` and path 
length `path`.
"""
function Eloss_µSi()
    Dminvec = [82.5, 317] * u"keV"
    DD01vec = [31.95, 66.4] * u"keV"
    pvec = [8.0, 7.0]
    dvec = [508, 1292] * u"keV"
    λvec = [131, 250] * u"keV"
    extrapolate(vec, k) = vec[1] * (1 - k) + vec[2] * k
    dpfunction = Δprob()
    p1GeV = MuscRat.muon_p(1u"GeV")

    x = [7.685074545223089e-5, 0.0004323905578247751, 0.001813403953654718, 0.006129759458666744, 0.016826296548409853, 0.03847041104177501, 0.07468542444861706, 0.1263986806126564, 0.1910928536006799, 0.2635700157200616, 0.33799263307399546, 0.40959903977735224, 0.4751623681658801, 0.5335239100304123, 0.5843968247098056, 0.6282086527148404, 0.665773347395409, 0.6980621439300688, 0.725688052810476, 0.7495258688857768, 0.7703381328859497, 0.7883654880543252, 0.8042391511197658, 0.8181875833892559, 0.8305565018526134, 0.8416043575016617, 0.8515305325734386, 0.8604439735792709, 0.8685034305438823, 0.875840624441936, 0.8825165876836277, 0.8886186589949316, 0.8942063692563273, 0.8993817240027229, 0.9041647037309595, 0.9086087968386626, 0.9127339807923409, 0.9165966500784148, 0.9201917076777463, 0.9235840178558824, 0.9267266240967387, 0.9297038802488505, 0.9324957361177313, 0.9351639847270156, 0.9377013790367099, 0.9400954416530567, 0.9423363876326569, 0.9444634688714055, 0.9464903809566951, 0.9484524636252607, 0.95028390260234, 0.9520253573732393, 0.9537068570171191, 0.9552920360297191, 0.9568099583293231, 0.9582804786539763, 0.9596966014812768, 0.961034581126346, 0.9623308618820874, 0.9635711617917846, 0.9647511580862334, 0.9658741477671121, 0.9669840782304073, 0.9680249398419165, 0.9690386814933922, 0.9700186637264927, 0.9709562182753173, 0.9718541966826775, 0.9727281400084127, 0.9735829301180585, 0.9743993350068829, 0.975189976009434, 0.9759552804025835, 0.9766993800445924, 0.9774232007020165, 0.9781358833726141, 0.9788059216901677, 0.9794594007561874, 0.9800862886551818, 0.9807054204183356, 0.9813010413473737, 0.9818798696892814, 0.9824307704875647, 0.9829742015164213, 0.9835053218198035, 0.9840285481072196, 0.9845305726738276, 0.9850244001912273, 0.9855007260237091, 0.9859686108652224, 0.9864092181694162, 0.9868494709247164, 0.9872737750404663, 0.9876948942768381, 0.988108546774116, 0.9885043172797564, 0.9888906786032188, 0.9892718595734014, 0.9896358040127534, 0.99]

    Enodes = LinRange(100, 508, 100)
    spl5 = interpolate(x, Enodes, FritschCarlsonMonotonicInterpolation())

    x = [0.00044905431092138483, 0.0018923459442988528, 0.00626060656646952, 0.01648807798251834, 0.03605130770871559, 0.06729060255544099, 0.11013660428109814, 0.162693882063288, 0.22141771904444651, 0.28262658594249757, 0.3432175155952201, 0.4007178968918101, 0.45404168758807983, 0.5024862657873438, 0.5459663054179461, 0.5849194712875303, 0.6194432400894483, 0.6500594148128065, 0.6772836337410655, 0.7015635744042535, 0.7231469013011615, 0.7425366404729608, 0.7599791180217552, 0.7755743500592152, 0.7897088842355917, 0.802560063162465, 0.8142572628293192, 0.8249470013093942, 0.8347343422375643, 0.8437101410314953, 0.8519505150885813, 0.8595921455189404, 0.866666305855122, 0.8732095661451684, 0.8793015726082283, 0.8850064581832107, 0.8903842658707705, 0.8953814719247568, 0.9000815741164448, 0.9045011756355285, 0.9086952164033705, 0.9126414906512395, 0.9164121248291256, 0.9199722620367459, 0.923309619327901, 0.9265013944133864, 0.9295939443038856, 0.9324537232729152, 0.9352219417670047, 0.9378887804999696, 0.9404272422429396, 0.9428519316815744, 0.9451514103222259, 0.947372723566025, 0.9494902287425595, 0.951499751064277, 0.9533939995520927, 0.9551907955892135, 0.9569306986559192, 0.958575261918474, 0.9601466083165308, 0.9616697234434864, 0.9631318192977886, 0.9645231866953196, 0.9658401553988252, 0.967104333696553, 0.9683230959017207, 0.9694711200540657, 0.9705844777119962, 0.9716582097369952, 0.97267482808239, 0.9736563935660246, 0.9745984106901973, 0.9755136821030654, 0.9763918162907294, 0.9772325344626338, 0.9780287284697327, 0.9787940242004101, 0.9795296034841935, 0.9802375804827929, 0.9809194728875457, 0.9815750276657201, 0.9822016518631306, 0.9828267018030039, 0.9834163544333663, 0.9839693976213686, 0.9845080073374765, 0.9850229683417082, 0.9855148387202295, 0.9859939120062017, 0.9864630336814869, 0.9869093450367831, 0.9873468446913252, 0.9877705795672608, 0.9881816875543015, 0.9885681276277627, 0.9889537348649455, 0.9893130792237423, 0.9896625376290481, 0.99]
    Enodes = LinRange(360, 1292, 100)
    spl15 = interpolate(x, Enodes, FritschCarlsonMonotonicInterpolation())


    function estimate_loss(p::Unitful.Momentum, path::Unitful.Length)
        k = ((path - 0.5u"mm") / 1u"mm") |> NoUnits
        loss_scale = (dpfunction(p, path) / dpfunction(p1GeV, path)) |> NoUnits 
        f = rand()
        loss_1GeV = 0.0u"keV"
        if f <= 0.01
            D0 = extrapolate(Dminvec, k)
            if D0 < 0u"keV"
                D0 *= 0.0
            end
            DeltaD01 = extrapolate(DD01vec, k)
            pwr = extrapolate(pvec, k)
            D1 = D0 + DeltaD01
            loss_1GeV = ((f / 0.01)^(1 / pwr) * (D1 - D0) + D0)
        elseif f >= 0.99
            c = 0.01
            d = extrapolate(dvec, k)
            λ = extrapolate(λvec, k)
            loss_1GeV = d - λ * log((1 - f) / c)
        else
            loss_1GeV = extrapolate([spl5(f), spl15(f)], k) * 1u"keV"
        end
        if loss_1GeV < 0.0u"keV"
            loss_1GeV *= 0
        end
        return loss_scale * loss_1GeV
    end
end