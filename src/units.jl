using Unitful

const mµ = 0.105659u"GeV/c^2"  # µ mass
const me = 0.510999u"MeV/c^2"  # e± mass

Ezero = 0u"GeV"
Pzero = 0u"GeV/c"
GeVc=1u"GeV/c"


@derived_dimension CRSpectrum (Unitful.𝐌/Unitful.𝐓^3)
@derived_dimension CRFlux inv(Unitful.𝐋^2*Unitful.𝐓)
Unitful.promote_unit(::S, ::T) where {S<:Unitful.EnergyUnits, T<:Unitful.EnergyUnits} = u"GeV"
Unitful.promote_unit(::S, ::T) where {S<:Unitful.MomentumUnits, T<:Unitful.MomentumUnits} = u"GeV/c"
