using Unitful

const me = 0.510999u"MeV/c^2"  # e± mass
const mµ = 105.659u"MeV/c^2"   # µ±mass
const mπ = 139.5706u"MeV/c^2"  # π± mass
const mp = 938.27208u"MeV/c^2" # proton mass
const mn =  939.5654133u"MeV/c^2" # neutron mass


Ezero = 0u"GeV"
Pzero = 0u"GeV/c"
GeVc=1u"GeV/c"


@derived_dimension CRSpectrum (Unitful.𝐌/Unitful.𝐓^3)
@derived_dimension CRFlux inv(Unitful.𝐋^2*Unitful.𝐓)
Unitful.promote_unit(::S, ::T) where {S<:Unitful.EnergyUnits, T<:Unitful.EnergyUnits} = u"GeV"
Unitful.promote_unit(::S, ::T) where {S<:Unitful.MomentumUnits, T<:Unitful.MomentumUnits} = u"GeV/c"
