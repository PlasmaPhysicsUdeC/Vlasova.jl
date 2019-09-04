"""
TODO
"""
function integrate!(plasma::Plasma, final_time::Real, dt::Real;
                    integrator::VlasovaIntegrator = integrator,
                    external_potential::Function = external_potential,
                    continue_from_backup::Bool = continue_from_backup,
                    save_distribution_times::Array{T} where T <: Real = save_distribution_times,
                    checkpoint_percent::Integer = checkpoint_percent,
                    velocity_filtering::Bool = velocity_filtering,
                    data_path::String = data_path,
                    FFTW_flags = FFTW_flags)

    println("Preparing integrator. This may take a while...")

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
    datasaver = DataSaver(plasma, Nt, dt, data_path, checkpoint_percent, continue_from_backup, save_distribution_times)

    # Do the magic!
    integrator(plasma, Nt, dt,
               poisson, external_potential,
               space_advection, velocity_advection,
               velocity_filtering,
               datasaver )

    return nothing;
end
