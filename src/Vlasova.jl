module Vlasova

import LinearAlgebra, FFTW, HDF5, Libdl, Dates


# Make compiled libraries available to julia
# push!(Base.DL_LOAD_PATH,"$(pwd())/src/fortran/sharedLibraries")

# Structs
export Specie, Plasma, Box
include("Box.jl")
include("Plasma.jl")

# Main integrator functions
export vlasova_integrator!
include("step.jl")
include("poisson.jl")
include("space_advection.jl")
include("velocity_advection.jl")
include("integrator.jl")

# Misc tools
include("tools/field_tools.jl")
include("tools/fourier_tools.jl")
include("tools/other_tools.jl")
include("tools/distributions.jl")

end
