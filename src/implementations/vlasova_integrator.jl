"""
        Takes a element of type Plasma and integrates it until final_time using steps of length dt.
        Accepts the optional ketwords:
            * continue_from_backup [false]: Tells whether the simulation should be restarted from a backup file.
            * checkpoint_percent [10]: Is the percentage accomplished between one flush of the data to the disk and another.
"""
function vlasova_integrator!(plasma, final_time, dt;
                             integrator::VlasovaIntegrator = verlet_velocity,
                             external_potential::Function = get_zero,
                             continue_from_backup::Bool = false,
                             save_distribution_times::Array{T} where T <: Real = Float64[],
                             checkpoint_percent::Integer = 10,
                             velocity_filtering::Bool = true,
                             progress_file::String = "/",
                             FFTW_flags = FFTW.ESTIMATE )         # TODO: Test the [nosave] case: checkpoint_percent = 100)
    # Number of iterations
    Nt = floor(Int, final_time/dt + 1)

    # Ensure saving distribution times are ordered Float64's
    save_distribution_times = sort( Float64.(save_distribution_times) )
    ## All saving times must be multiple of dt
    @assert all( isinteger.( round.(save_distribution_times ./ dt, digits = 10) )
                 ) "Not all times to save the distribution function are multiples of dt"
    @assert isapprox(save_distribution_times[end] - final_time, 0,
                     atol = 1e-10) "The times to save the distribution function can not be larger than final_time"
    ## Always save the last instant
    (final_time in save_distribution_times) ? nothing : append!(save_distribution_times, final_time)
    

    
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
    ## Also, initialize h5 files (saving first checkpoint) or restore data and allow to save the whole DF at save_distribution times
    datasaver = DataSaver(plasma, Nt, checkpoint_percent, continue_from_backup, save_distribution_times)
        
    # Do the magic!
    notify("Everything initialized.")
    integrator(plasma, Nt, dt,
               poisson, external_potential,
               space_advection, velocity_advection,
               velocity_filtering,
               datasaver,
               progress_file = progress_file)
    return; # nothing
end
