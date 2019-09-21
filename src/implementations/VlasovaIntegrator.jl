"""
```julia
redefine_integrator( codeblock::Expr )
```

Inject `codeblock` explicitly into the main integrator function, such
that it gets executed into every interation of the main temporal loop
(once for every time step).

# Notes
* You should only use this function if you **really** know what you are doing.
  This function is very likely to break the code.
"""
redefine_integrator( codeblock::Expr ) = @eval begin
    function (integrator::VlasovaIntegrator)(plasma::Plasma,
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
        progressbar = ProgressMeter.Progress(length(iteration_axis),
                                             barglyphs = ProgressMeter.BarGlyphs("[=> ]"),
                                             output = buff )

        # Sync progressbar when continuing from a backup
        progressbar.counter = datasaver.last_iteration_saved

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
                            #get_gradient_correction!(grad, poisson!, electricfield) # Get $ grad = \nabla |E|^2 $
                            grad[1] = (2 * electricfield[1] .* (chargedensity .+ 1)) # TODO: This line works (1d), but idk why
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

            # Inject external code
            $codeblock
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
        return nothing;
    end
end

# Define the integrator function injecting an empty code block
redefine_integrator(quote end)
