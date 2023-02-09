module MuscRat

export Line, volume, Box, HCylinder
export path, smallest_radius
export µspectrum, µspectrum_p, µspectrum_reyna, µspectrum_reyna_p
export Eloss_function

include("geometry.jl")
include("cr_spectrum.jl")
include("energy_loss.jl")

greet() = print("Hello World!")

end # module MuscRat
