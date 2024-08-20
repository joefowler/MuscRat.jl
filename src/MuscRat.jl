module MuscRat

export Line, Box, HCylinder, VCylinder, Cylinder, Sphere
export path, volume, smallest_radius
export µspectrum, µspectrum_p, µspectrum_reyna, µspectrum_reyna_p
export CRGenerator, CRMuonGenerator, ParmaGenerator, generate
export Eloss_functions, dEdx, Δprob

include("geometry.jl")
include("units.jl")
include("parma.jl")
include("cr_spectrum.jl")
include("cr_generator.jl")
include("energy_loss.jl")

end # module MuscRat
