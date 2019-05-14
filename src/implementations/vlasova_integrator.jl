"""
        Takes a element of type Plasma and integrates it until final_time using steps of length dt.
        Accepts the optional ketwords:
            * continue_from_backup [false]: Tells whether the simulation should be restarted from a backup file.
            * checkpoint_percent [10]: Is the percentage accomplished between one flush of the data to the disk and another.
"""
function vlasova_integrator!(plasma, final_time, dt;
                             integrator::VlasovaIntegrator = verlet_velocity,
                             continue_from_backup::Bool = false,
                             checkpoint_percent::Integer = 10,
                             velocity_filtering::Bool = true,
                             progress_file::String = "/",
                             FFTW_flags = FFTW.ESTIMATE )         # TODO: Test the [nosave] case: checkpoint_percent = 100)
    # Number of iterations
    Nt = floor(Int, final_time/dt + 1)
    
    # Make fortran libraries available to Julia
    push!(Libdl.DL_LOAD_PATH, joinpath(dirname(@__FILE__), "../../deps/usr/lib"))     # TODO: This is not the best way to accomplish that
    
    # Initialize objects
    notify("Preparing integrator data. This may take a while...")
    ##  To solve poisson equation
    poisson = Poisson(plasma, FFTW_flags = FFTW_flags)
    ## Advections
    space_advection = SpaceAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
    velocity_advection = VelocityAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
    ## To save data in memory and flush it to disk on checkpoints
    ## Also, initialize h5 files (saving first checkpoint) or restore data
    datasaver = DataSaver(plasma, Nt, checkpoint_percent, continue_from_backup)
        
    # Do the magic!
    notify("Everything initialized.")
    integrator(plasma, Nt, dt,
               poisson, space_advection, velocity_advection,
               velocity_filtering,
               datasaver,
               progress_file = progress_file)
    return; # nothing
end
