function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                         time_manager::TimeManager,
                                         poisson!::Poisson,
                                         external_potential::Function,
                                         space_advection::SpaceAdvection,
                                         velocity_advection::VelocityAdvection,
                                         velocity_filtering::Bool,
                                         datasaver::DataSaver;
                                         progress_file::String)

    # Preallocate to make operations in place
    chargedensity = get_density( plasma )
    electricfield = poisson!( chargedensity )
    prop = [ similar(electricfield[d], Complex{Float64})
             for d in axes(electricfield, 1) ]
    if 'C' in integrator.sequence
        grad = deepcopy( electricfield )
    else
        grad = nothing
    end

    # Iteration axis
    iteration_axis = (datasaver.last_iteration_saved + 1):time_manager.final_iteration    
    notify("Entering main loop... $(Dates.now())", filename = progress_file, mode = "w")

    # Start counting time
    time_manager( start = Dates.now() )

    # Go!
    for t in iteration_axis
        time = (t-2) * time_manager.dt
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
                    get_gradient_correction!(grad, poisson!, electricfield) # Get $ grad = \nabla |E|^2 $
                end

                velocity_advection(plasma, electricfield, grad, prop,
                                   advection_number = vel_adv_num,
                                   gradient_number = grad_adv_num,
                                   is_gradient_advection = isC,
                                   filtering = velocity_filtering && (vel_adv_num == 1 )) # Apply filter just once per time iteration
            end
        end

        # Save data and notify progress at checkpoints
        datasaver(plasma, time_manager, t, progress_file = progress_file)
    end
    return 0;
end
