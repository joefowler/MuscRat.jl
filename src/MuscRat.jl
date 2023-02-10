module MuscRat

export Line, Box, HCylinder, Sphere
export path, volume, smallest_radius
export µspectrum, µspectrum_p, µspectrum_reyna, µspectrum_reyna_p
export Eloss_function
export CRMuonGenerator, generate

include("geometry.jl")
include("cr_generator.jl")
include("cr_spectrum.jl")
include("energy_loss.jl")

end # module MuscRat
