using Vlasova

function main(box::Box,                                    # Box parameters
              species::Array{Specie},                      # Initial conditions
              dt, final_time;                              # Final time
              continue_from_backup = continue_from_backup, # Backup?
              FFTW_flags = Vlasova.FFTW.ESTIMATE )         # FFTW
    
    # Integrate until
    Nt = floor(Int, final_time/dt + 1)

    # Create whole plasma
    plasma = Plasma( species, box )

    vlasova_integrator!(plasma, Nt, dt,
                        continue_from_backup = continue_from_backup,
                        FFTW_flags = FFTW_flags )
end

# Get parameters
include("parameters.jl")

# Continue from backup ? 
(@isdefined continue_from_backup) ? nothing : continue_from_backup = false

# FFTW flags
(@isdefined FFTW_flags) ? nothing : FFTW_flags = Vlasova.FFTW.ESTIMATE

# Perform the simulation
simulation_path = "data/$simulation_name"
mkpath(simulation_path)
run(`cp parameters.jl $simulation_path`)
Vlasova.notify("Folder $simulation_path prepared")

# Allow multithreading
vlasova_multithread(FFTW_NUM_THREADS = num_threads,
                    OMP_NUM_THREADS = num_threads  )

# Run!
main(box,
     species,
     dt, final_time,
     continue_from_backup = continue_from_backup,
     FFTW_flags = FFTW_flags )
