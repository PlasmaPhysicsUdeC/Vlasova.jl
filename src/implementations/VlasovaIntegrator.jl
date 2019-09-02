function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                         Nt::Integer, dt::Float64,
                                         poisson!::Poisson,
                                         external_potential::Function,
                                         space_advection::SpaceAdvection,
                                         velocity_advection::VelocityAdvection,
                                         velocity_filtering::Bool,
                                         datasaver::DataSaver,
                                         outputs::Array{T, 1}) where T <: IO

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
    progressbars = [ ProgressMeter.Progress(length(iteration_axis),
                                            output = op ) for op in outputs ]

    # Sync progressbars when continuing from a backup
    [ProgressMeter.next!.(progressbars) for i in 2:datasaver.last_iteration_saved ]

    # Go!
    println.(outputs, "Starting integration @ $(Dates.now())"); flush.(outputs)
    for t in iteration_axis
        time = (t-2) * dt
        pos_adv_num = 0
        vel_adv_num = 0
        grad_adv_num = 0
        for a in integrator.sequence
            if a == 'A'
                pos_adv_num += 1
                space_advection(plasma, advection_number = pos_adv_num)
                get_density!(chargedensity, plasma)

                time += space_advection.coefficients[ pos_adv_num ] # Updtate time after advection
                poisson!(electricfield, chargedensity, external_potential = external_potential( time, plasma.box ) )

            else # Velocity advection (B or C)
                vel_adv_num += 1

                isC = ( a == 'C' )
                if isC
                    grad_adv_num += 1
                    #get_gradient_correction!(grad, poisson!, electricfield) # Get $ grad = \nabla |E|^2 $
                    grad[1] = (2 * electricfield[1] .* (chargedensity .+ 1)) # TODO: This line works (1d), but idk why
                end

                velocity_advection(plasma, electricfield, grad, prop,
                                   advection_number = vel_adv_num,
                                   gradient_number = grad_adv_num,
                                   is_gradient_advection = isC,
                                   filtering = velocity_filtering ) # Apply filter at all velocity advections
            end
        end

        # Save data and update progressbar
        datasaver(plasma, t)
        ProgressMeter.next!.(progressbars)
    end
    return 0;
end
