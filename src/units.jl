using Unitful

const mยต = 0.105659u"GeV/c^2"  # ยต mass
const me = 0.510999u"MeV/c^2"  # eยฑ mass

Ezero = 0u"GeV"
Pzero = 0u"GeV/c"
GeVc=1u"GeV/c"


@derived_dimension CRSpectrum (Unitful.๐/Unitful.๐^3)
@derived_dimension CRFlux inv(Unitful.๐^2*Unitful.๐)
Unitful.promote_unit(::S, ::T) where {S<:Unitful.EnergyUnits, T<:Unitful.EnergyUnits} = u"GeV"
Unitful.promote_unit(::S, ::T) where {S<:Unitful.MomentumUnits, T<:Unitful.MomentumUnits} = u"GeV/c"
