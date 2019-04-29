"""
        Takes a element of type Plasma and integrates it over Nt time steps of length dt.
        Accepts the optional ketwords:
            * continue_from_backup (false): Tells whether the simulation should be restarted from a backup file.
            * checkpoint_percent (10): Is the percentage accomplished between one flush of the data to the disk and another.
        """
function vlasova_integrator!(plasma, Nt, dt;
                             continue_from_backup::Bool = false,
                             checkpoint_percent::Integer = 5,
                             velocity_filtering::Bool = true,
                             FFTW_flags = FFTW.ESTIMATE )         # TODO: Test the [nosave] case: checkpoint_percent = 100)
    
    # checkpoints
    @assert 1 <= checkpoint_percent <= 100  "checkpoint_percent must be valued between 1 and 100"
    @assert isinteger(100/checkpoint_percent) "checkpoint_percent must be a whole divisor of 100"
    checkpoint_axis = collect(Int, range(1, stop = Nt, length = Int(100/checkpoint_percent) + 1 ))
    checkpoint_step = checkpoint_axis[2] - checkpoint_axis[1]

    # Make fortran libraries available to Julia
    push!(Libdl.DL_LOAD_PATH, joinpath(dirname(@__FILE__), "../deps/usr/lib"))     # TODO: Try to change this
    
    # Initialize objects
    poisson! = Poisson(plasma, FFTW_flags = FFTW_flags)
    space_advection! = SpaceAdvection(plasma, dt, FFTW_flags = FFTW_flags)
    velocity_advection! = VelocityAdvection(plasma, dt, FFTW_flags = FFTW_flags)
    
    # Initialize quantities
    chargedensity = Array{Float64}(undef, plasma.box.Nx)
    electricfield = poisson!(chargedensity)
    
    ## Timed quantities
    timed_chargedensity = Array{Float64}(undef, (plasma.box.Nx..., checkpoint_step + 1) )
    timed_kinen = Array{Float64}(undef, (checkpoint_step + 1, plasma.number_of_species) )

    last_iteration_saved = Ref(0)
    if continue_from_backup
        restore_variables!( plasma, last_iteration_saved )
        notify("Integrator variables restored from backup")
        
        # Starting conditions            
        get_density!(chargedensity, plasma )
        poisson!(electricfield, chargedensity)
    else
        initialize_h5files( plasma, Nt)
        notify("Integrator variables and files initialized")

        # Starting conditions
        get_density!(chargedensity, plasma )
        poisson!(electricfield, chargedensity)
        timed_chargedensity[plasma.box.space_axis, 1] .= chargedensity
        timed_kinen[1, :] .= get_kinetic_energies( plasma )
        flushdata(plasma, timed_chargedensity, timed_kinen, last_iteration_saved, 1, Nt)

        # First velocity advection
        velocity_advection!( plasma, electricfield )
    end
    iteration_axis = (last_iteration_saved[] + 1):Nt # Is this correct when restoring variables?
    
    # Main temporal loop
    notify("\tEntering main loop...")
    start_time = Dates.now() # Measure execution time of the main loop
    for t in iteration_axis[2:end]
        step!( plasma, chargedensity, electricfield,
               poisson!, space_advection!, velocity_advection!,
               stepnumber = t, velocity_filtering = velocity_filtering)
        
        # Save every instant
        inst = t - last_iteration_saved[]
        timed_chargedensity[plasma.box.space_axis, inst] .= chargedensity
        timed_kinen[inst, :] .= get_kinetic_energies( plasma )
        if t in checkpoint_axis
            flushdata( plasma, timed_chargedensity, timed_kinen, last_iteration_saved, t, Nt)
            last_checkpoint_saved = findfirst( t .== checkpoint_axis ) - 1 # TODO: use an accumulator for this
            
            elapsed_time = round(Dates.now() - start_time, Dates.Second )
            notify("\t $(last_checkpoint_saved*checkpoint_percent)% accomplished in $elapsed_time")
        end
    end
    return 0;
end

function restore_variables!( plasma, last_iteration_saved )
    # Read data from H5 files to continue the simulation from the last checkpoint
    for s in 1:plasma.number_of_species
        species_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/"*plasma.species[s].name*".h5", "r")
        s == 1 ? (last_iteration_saved[] = read( species_file["last_iteration_saved"])[1]) : nothing
        plasma.species[s].distribution .= read( species_file["distribution"] )
        
        HDF5.close(species_file)
    end
    return 0;
end

function initialize_h5files(plasma, Nt)
    # File for the charge density
    shared_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/shared_data.h5", "w")
    shared_file["chargedensity"] = Array{Float64}(undef, (plasma.box.Nx..., Nt))
    shared_file["total_kinetic_energy"] = Array{Float64}(undef, Nt)
    shared_file["specie_names"] = Array([plasma.species[s].name for s in plasma.specie_axis])
    HDF5.close(shared_file)

    # File for the distribution functions
    for s in 1:plasma.number_of_species
        species_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/"*plasma.species[s].name*".h5", "w")
        s == 1 ? ( species_file["last_iteration_saved"] = [0] ) : nothing
        species_file["distribution"] = Array{Float64}(undef, plasma.box.N )
        HDF5.close(species_file)
    end
end

function flushdata(plasma, timed_chargedensity, timed_kinen, last_iteration_saved, t, Nt)
    #=
    Flushed saved data to files
    =#

    last_it = last_iteration_saved[]
    #TODO: turn space_axis into UnitRanges
    
    shared_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/shared_data.h5", "r+")
    shared_file["chargedensity"][UnitRange.(1, plasma.box.Nx)..., (last_it+1):t] = timed_chargedensity[plasma.box.space_axis, 1:(t-last_it)]
    shared_file["total_kinetic_energy"][(last_it+1):t] = dropdims( sum( timed_kinen[1:(t-last_it), :], dims = 2 ), dims = 2 )
    
    HDF5.close(shared_file)

    # File for the distribution functions
    for s in plasma.specie_axis
        species_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/"*plasma.species[s].name*".h5", "r+")
        s == 1 ? (species_file["last_iteration_saved"][:] = t) : nothing
        species_file["distribution"][UnitRange.(1, plasma.box.N)...] = plasma.species[s].distribution
        HDF5.close(species_file)
    end

    # Update saving indices
    last_iteration_saved[] = t

    if t == Nt  # At last iteration, specify that the files are now complete.
        shared_file = HDF5.h5open("data/"*plasma.box.simulation_name*"/shared_data.h5", "r+")
        HDF5.attrs(shared_file)["status"] = "Complete"
        HDF5.close(shared_file)
    end
    return 0;
end
