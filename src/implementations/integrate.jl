"""
```julia
integrate!(plasma::Plasma, final_time::Real, dt::Real; kwargs...)
```

Evolve `plasma` in place from `t = 0` up to `final_time` with time step `dt`, solving the
Vlasov-Poisson system with a symplectic split-step scheme.

# Keyword arguments
Every keyword defaults to the Vlasova variable of the same name (see [Default behavior](@ref)):
* `integrator`: [`VlasovaIntegrator`](@ref) used to split each time step.
* `external_potential`: function adding an external electrostatic potential.
* `continue_from_backup`: resume an unfinished simulation stored in `data_path`.
* `save_distribution_times`: times at which the full distribution function is saved.
* `checkpoint_percent`: percentage of the simulation between checkpoints saved to disk.
* `velocity_filtering`: apply anisotropic filtering in the transformed-velocity space.
* `data_path`: folder where results are saved. If empty, nothing is written to disk.
* `FFTW_flags`: flags passed to the FFTW plans.
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
