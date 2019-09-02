"""
        Take an element of type Plasma and integrate it until final_time using steps of length dt.
        Accepts the optional keywords:
            * continue_from_backup [false]: Tells whether the simulation should be restarted from a backup file.
            * checkpoint_percent [10]: Is the percentage accomplished between one flush of the data to the disk and another.
"""
function integrate(plasma, final_time, dt;
                   integrator::VlasovaIntegrator = integrator,
                   external_potential::Function = external_potential,
                   continue_from_backup::Bool = continue_from_backup,
                   save_distribution_times::Array{T} where T <: Real = save_distribution_times,
                   checkpoint_percent::Integer = checkpoint_percent,
                   velocity_filtering::Bool = velocity_filtering,
                   save_path::String = save_path,
                   FFTW_flags = FFTW_flags)

    if (size(save_distribution_times, 1) !== 0)
        @assert (save_path !== "/") "The variable save_path must be specified if you want to save distributions"
    end
    if ( checkpoint_percent < 100 )
        @assert (save_path !== "/") "The variable save_path must be specified if you want to use checkpoints"
    end
    mkpath(save_path)
    save_data = (save_path !== "/") ||
        (size(save_distribution_times, 1) !== 0) ||
        ( checkpoint_percent < 100 )

    # Output streams
    outputs = [ stdout ]
    if save_data
        fid = open( joinpath(save_path, "progress_file"), "w")
        outputs = [stdout, fid]
    end
    println.(outputs, "Preparing integrator. This may take a while..."); flush.(outputs)

    # Number of time iterations
    Nt = round(Int, final_time / dt) + 1

    # Multithread
    (NUM_THREADS > 1) ? FFTW.set_num_threads( NUM_THREADS ) : nothing

    # Initialize objects
    ##  To solve poisson equation
    poisson = Poisson(plasma, FFTW_flags = FFTW_flags)
    ## Advections
    space_advection = SpaceAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
    velocity_advection = VelocityAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
    ## To save data in memory and flush it to disk on checkpoints
    ## Also, initialize h5 files (saving first checkpoint) or restore data and allow to save the whole DF at save_distribution times
    datasaver = DataSaver(plasma, Nt, dt, save_data, save_path, checkpoint_percent, continue_from_backup, save_distribution_times)

    # Do the magic!
    integrator(plasma, Nt, dt,
               poisson, external_potential,
               space_advection, velocity_advection,
               velocity_filtering,
               datasaver,
               outputs)

    # close progress_file if it corresponds
    save_data ? close(fid) : nothing
    return nothing;
end
