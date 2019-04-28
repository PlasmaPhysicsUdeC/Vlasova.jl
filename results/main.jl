using Vlasova

function main(box::Box,     # Box
              species::Array{Specie},
              dt, final_time;                              # Final time
              continue_from_backup = continue_from_backup, # Backup?
              FFTW_flags = Vlasova.FFTW.ESTIMATE)          # FFTW
    
    # Integrate until
    Nt = floor(Int, final_time/dt + 1)

    # Create whole plasma
    plasma = Plasma( species, box )

    vlasova_integrator!(plasma, Nt, dt, continue_from_backup = continue_from_backup)
end

# take parameters
include("parameters.jl")

# Continue from backup ? 
(@isdefined continue_from_backup) ? nothing : continue_from_backup = false

# FFTW flags
(@isdefined FFTW_flags) ? nothing : FFTW_flags = Vlasova.FFTW.ESTIMATE

# Multithreading
## For FFTW
(@isdefined FFTW_NUM_THREADS) ? Vlasova.FFTW.set_num_threads(FFTW_NUM_THREADS) : nothing
## For advections (FORTRAN)
(@isdefined OMP_NUM_THREADS) ? (ENV["OMP_NUM_THREADS"] = OMP_NUM_THREADS) : nothing 

# Perform the simulation
simulation_path = "data/$simulation_name"
mkpath(simulation_path)
run(`cp parameters.jl $simulation_path`)
Vlasova.notify("Folder $simulation_path prepared")

main(box,
     species,
     dt, final_time,
     continue_from_backup = continue_from_backup,
     FFTW_flags = FFTW_flags)
