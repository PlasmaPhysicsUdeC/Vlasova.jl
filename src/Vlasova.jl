module Vlasova

import Statistics, LinearAlgebra, FFTW, HDF5, Dates, QuadGK, ProgressMeter

# Structs
export Box, VlasovaIntegrator, Specie, Plasma, get_zero
include("structs/01_Box.jl")
include("structs/01_Zero.jl")
include("structs/01_VlasovaIntegrator.jl")
include("structs/01_Specie.jl")
include("structs/02_Plasma.jl")
## Private
include("structs/03_DataSaver.jl")
include("structs/03_Poisson.jl")
include("structs/03_SpaceAdvection.jl")
include("structs/03_VelocityAdvection.jl")

# Function implementations
include("implementations/Plasma.jl")
include("implementations/Poisson.jl")
include("implementations/DataSaver.jl")
include("implementations/SpaceAdvection.jl")
include("implementations/VelocityAdvection.jl")
include("implementations/VlasovaIntegrator.jl")

# Main integrator
export integrate
include("implementations/integrate.jl")

# Misc tools
# Multiple exports here. They are specified in each file.
include("tools/integrators.jl")
include("tools/field_tools.jl")
include("tools/fourier_tools.jl")
include("tools/other_tools.jl")
include("tools/distributions.jl")
include("tools/defaults.jl")

end
