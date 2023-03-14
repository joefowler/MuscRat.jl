module MuscRat

export Line, Box, HCylinder, VCylinder, Cylinder, Sphere
export path, volume, smallest_radius
export µspectrum, µspectrum_p, µspectrum_reyna, µspectrum_reyna_p
export CRGenerator, CRMuonGenerator, generate
export Eloss_function

include("geometry.jl")
include("units.jl")
include("cr_spectrum.jl")
include("cr_generator.jl")
include("energy_loss.jl")

end # module MuscRat
