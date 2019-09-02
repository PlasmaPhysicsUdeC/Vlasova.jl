using Vlasova, Statistics

# Get parameters
parametersfile = "parameters.jl"

is_from_backup = (@isdefined continue_from_backup) && continue_from_backup
if !isempty( ARGS ) && !is_from_backup
    parametersfile = ARGS[1]
end

println("Parameters file: ", parametersfile)
include( parametersfile )

Base.Sys.set_process_title("Vlasova-" )    # Set process title of the simulation #TODO: append name

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
integrate(plasma, final_time, dt )
