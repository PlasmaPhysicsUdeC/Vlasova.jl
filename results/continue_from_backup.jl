"""
    This file is a parti of the Vlasova package, developed by Jorge Gidi for 
    the group of plasma physics of the University of Concepción
"""

import HDF5

# Detect if the argument given is a string
@assert (!isempty(ARGS) && length(ARGS)==1 ) "You must start the program passing the simulation name as the first (and only) argument.
                                               Example: user@computer \$ julia continue_from_backup simulation_name"

simulation_path = "data/"*ARGS[1]

@assert isfile(simulation_path*"/shared_file.h5") "There is no file shared_file.h5 into $simulation_path/. That file should contain data previously calculated by the simulation."

# Detect if file is already completed
file = HDF5.h5open(simulation_path*"/shared_file.h5", "r")
iscomplete = HDF5.exists(HDF5.attrs(file), "status")
HDF5.close(file)

@assert !is_complete "This simulation has already ended successfully"

# Restore parameters.jl to main directory
run(`cp $simulationPath/parameters.jl parameters.jl`)

# Re-start simulation
const continue_from_backup = true
include("main.jl")
