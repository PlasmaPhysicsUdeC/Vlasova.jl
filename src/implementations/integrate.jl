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

    # Multithread FFTW plans
    # Threading for non-FFT operations is managed through Threads.@threads and Strided.@strided
    FFTW.set_num_threads( Threads.nthreads() )

    # Initialize objects
    TimerOutputs.@timeit_debug timer "Initialize objects" begin
        ##  To solve poisson equation
        TimerOutputs.@timeit_debug timer "Poisson" poisson = Poisson(plasma, FFTW_flags = FFTW_flags)
        ## Advections
        TimerOutputs.@timeit_debug timer "SpaceAdvection" space_advection = SpaceAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
        TimerOutputs.@timeit_debug timer "VelocityAdvection" velocity_advection = VelocityAdvection(plasma, integrator, dt, FFTW_flags = FFTW_flags)
        ## To save data in memory and flush it to disk on checkpoints
        ## Also, initialize h5 files (saving first checkpoint) or restore data and allow to save the whole DF at save_distribution times
        TimerOutputs.@timeit_debug timer "DataSaver" datasaver = DataSaver(plasma, Nt, dt, data_path, checkpoint_percent, continue_from_backup, save_distribution_times)
    end
    # Do the magic!
    TimerOutputs.@timeit_debug timer "Main integrator" begin
        integrator(plasma, Nt, dt,
                   poisson, external_potential,
                   space_advection, velocity_advection,
                   velocity_filtering,
                   datasaver )
    end

    return nothing;
end
