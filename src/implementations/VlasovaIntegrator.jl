"""
```julia
inject_to_integrator(;before_loop::Expr, inside_loop::Expr, after_loop::Expr )
```

Inject blocks of code explicitly into the main integrator function.
This function requires three keyword arguments:
* `before_loop`: Executed just before entering the integration loop.
* `inside_loop`: Executed at the end of every iteration of the integration loop, but
   before saving the data. This block will be executed after every `dt` is accomplished.
* `after loop`: Executed just after ending the integration loop. It is the last thing to do
   when the integrator is called.

# Notes
* You should only use this function if you **really** know what you are doing.
  This function is very likely to break the code.
"""
function inject_to_integrator(;before_loop::Expr, inside_loop::Expr, after_loop::Expr )
    # If the injected code is not empty, debug it with TimerOutputs
    if length(before_loop.args) > 1
        before_loop = quote
            TimerOutputs.@timeit_debug timer "[injected] before loop" $before_loop
        end
    end
    if length(inside_loop.args) > 1
        inside_loop = quote
            TimerOutputs.@timeit_debug timer "[injected] inside loop" $inside_loop
        end
    end
    if length(after_loop.args) > 1
        after_loop = quote
            TimerOutputs.@timeit_debug timer "[injected] after loop" $after_loop
        end
    end

    # Define integrator function
    @eval function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                                   Nt::Integer, dt::Float64,
                                                   poisson!::Poisson,
                                                   external_potential::Function,
                                                   space_advection::SpaceAdvection,
                                                   velocity_advection::VelocityAdvection,
                                                   velocity_filtering::Bool,
                                                   datasaver::DataSaver )

        # Preallocate to make operations in place
        chargedensity = get_density( plasma )
        electricfield = poisson!( chargedensity )
        prop = [ similar(electricfield[d], Complex{Float64})
                 for d in axes(electricfield, 1) ]
        if 'C' in integrator.sequence
            grad = [ similar( electricfield[d] )
                     for d in axes(electricfield, 1) ]
        else
            grad = nothing
        end

        # Iteration axis
        iteration_axis = (datasaver.last_iteration_saved + 1):Nt

        # Make progress indicators
        buff = IOBuffer()
        progressbar = ProgressMeter.Progress(Nt,
                                             barglyphs = ProgressMeter.BarGlyphs("[=> ]"),
                                             output = buff )

        # Sync progressbar when continuing from a backup
        progressbar.counter = datasaver.last_iteration_saved

        # Inject code before the main integration loop.
        $before_loop

        # Go!
        start_str = "Starting integration @ $(Dates.now())\n"
        print( start_str )
        for t in iteration_axis
            time = (t-2) * dt
            pos_adv_num = 0
            vel_adv_num = 0
            grad_adv_num = 0
            for a in integrator.sequence
                if a == 'A'
                    pos_adv_num += 1
                    TimerOutputs.@timeit_debug timer "space advection" begin
                        space_advection(plasma, advection_number = pos_adv_num)
                    end

                    TimerOutputs.@timeit_debug timer "get density" begin
                        get_density!(chargedensity, plasma)
                    end

                    time += space_advection.coefficients[ pos_adv_num ] # Updtate time after advection
                    TimerOutputs.@timeit_debug timer "get electric field" begin
                        poisson!(electricfield,
                                 chargedensity,
                                 external_potential = external_potential( time, plasma.box ) )
                    end

                else # Velocity advection (B or C)
                    vel_adv_num += 1

                    isC = ( a == 'C' )
                    if isC
                        grad_adv_num += 1
                        TimerOutputs.@timeit_debug timer "get force gradient" begin
                            get_gradient_correction!(grad, poisson!, electricfield) # Get $ grad = \nabla |E|^2 $
                            #grad[1] = (2 * electricfield[1] .* (chargedensity .+ 1)) # TODO: This line works (1d), but idk why
                        end
                    end

                    TimerOutputs.@timeit_debug timer "velocity advection" begin
                        velocity_advection(plasma, electricfield, grad, prop,
                                           advection_number = vel_adv_num,
                                           gradient_number = grad_adv_num,
                                           is_gradient_advection = isC,
                                           filtering = velocity_filtering ) # Apply filter at all velocity advections
                    end
                end
            end

            # Inject code at the end of every iteration of the integration loop.
            $inside_loop

            # Save data
            TimerOutputs.@timeit_debug timer "save data" datasaver(plasma, t)
            # Update progressbar in buffer
            ProgressMeter.next!(progressbar)
            buff_str = String(take!(buff))
            ## Buffer to stdout
            print( buff_str )
            ## Buffer to progressfile only if progressbar changed
            if datasaver.save_data && !isempty( buff_str )
                write( joinpath(datasaver.path, "progressfile"),
                       start_str * buff_str[2:end-3] )
            end
        end

        # Inject code at the end of the integration loop.
        $after_loop

        return nothing;
    end
end

# Define the integrator function injecting empty code blocks.
inject_to_integrator(before_loop = quote end,
                     inside_loop = quote end,
                     after_loop = quote end )
