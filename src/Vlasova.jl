module Vlasova

import LinearAlgebra, FFTW, HDF5, Libdl, Dates

# Structs
export Box, VlasovaIntegrator, Specie, Plasma
include("structs/01_Box.jl")
include("structs/01_VlasovaIntegrator.jl")
include("structs/01_Specie.jl")
include("structs/02_Plasma.jl")
## Private
include("structs/03_DataSaver.jl")
include("structs/03_Poisson.jl")
include("structs/03_SpaceAdvection.jl")
include("structs/03_VelocityAdvection.jl")

# Function implementations
include("implementations/VlasovaIntegrator.jl")
include("implementations/Plasma.jl")
include("implementations/DataSaver.jl")
include("implementations/Poisson.jl")
include("implementations/SpaceAdvection.jl")
include("implementations/VelocityAdvection.jl")

# Main integrator
export vlasova_integrator!
include("implementations/vlasova_integrator.jl")

# Extras
export vlasova_multithread
include("extras/multithreading.jl")

# Misc tools
include("tools/integrators.jl")
include("tools/field_tools.jl")
include("tools/fourier_tools.jl")
include("tools/other_tools.jl")
include("tools/distributions.jl")

end
