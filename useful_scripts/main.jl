using Vlasova

# Get parameters
parametersfile = "parameters.jl"

is_from_backup = (@isdefined continue_from_backup) && continue_from_backup
if !isempty( ARGS ) && !is_from_backup
    parametersfile = ARGS[1]
end

println("Parameters file: ", parametersfile)
include( parametersfile )

Base.Sys.set_process_title("julia "*simulation_name )    # Set process title of the simulation

# Continue from backup ? 
(@isdefined continue_from_backup) ? nothing : continue_from_backup = false
# FFTW flags
(@isdefined FFTW_flags) ? nothing : FFTW_flags = Vlasova.FFTW.ESTIMATE
# Multithreading
(@isdefined num_threads) ? vlasova_multithread( num_threads  ) : nothing
# Integrator
(@isdefined integrator) ? nothing : integrator = verlet_velocity
# Velocity filtering
(@isdefined velocity_filtering) ? nothing : velocity_filtering = true
# External potential
(@isdefined external_potential) ? nothing : external_potential = Vlasova.get_zero
# Save distribution
(@isdefined save_distribution_times) ? nothing : save_distribution_times = Float64[]
# Checkpoint percent
(@isdefined checkpoint_percent) ? nothing : checkpoint_percent = 10

# Prepare the simulation folder
simulation_path = "data/$simulation_name"           # Set name to save data
mkpath(simulation_path)                             # Create folder
run(`cp $parametersfile $simulation_path/parameters.jl`)            # Copy parameters.jl to simulation folder
println("Folder $simulation_path prepared")
progress_file = simulation_path*"/progress_file"

# Run!
# Create Specie array
species = [ Specie(name = name[s],
                   charge = charge[s],
                   mass = mass[s],
                   temperature = temperature[s],
                   distribution = initial_distribution(box, perturbate = perturbed[s]) ) for s in 1:length(name) ]
# Create Plasma
plasma = Plasma(species, box)

# Go!
integrate(plasma, final_time, dt,
                    save_distribution_times = save_distribution_times,
                    continue_from_backup = continue_from_backup,
                    external_potential = external_potential,
                    velocity_filtering = velocity_filtering,
                    checkpoint_percent = checkpoint_percent,
                    integrator = integrator,
                    FFTW_flags = FFTW_flags,
                    progress_file = progress_file)
