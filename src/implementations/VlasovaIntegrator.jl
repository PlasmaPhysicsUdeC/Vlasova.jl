
function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                         time_manager::TimeManager,
                                         poisson!::Poisson,
                                         external_potential::Function,
                                         sadv!::SpaceAdvection,
                                         vadv!::VelocityAdvection,
                                         velocity_filtering::Bool,
                                         datasaver::DataSaver;
                                         progress_file::String)

    # Preallocate to make operations in place
    chargedensity = get_density( plasma )
    electricfield = poisson!( chargedensity )

    # Iteration axis
    iteration_axis = (datasaver.last_iteration_saved + 1):time_manager.final_iteration
    
    notify("Entering main loop...", filename = progress_file, mode = "w")

    # Start counting time
    time_manager( start = Dates.now() )

    # Go!
    for t in iteration_axis
        time = (t-2) * time_manager.dt
        pos_adv_num = 0
        vel_adv_num = 0
        # TODO: Allow for merging
        for a in integrator.sequence
            if a == 'A'
                pos_adv_num += 1
                sadv!(plasma, advection_number = pos_adv_num)
                get_density!(chargedensity, plasma)
                
                time += sadv!.coefficients[ pos_adv_num ] # Updtate time after advection
                poisson!(electricfield, chargedensity, external_potential = external_potential( time, plasma.box ) )
            else
                vel_adv_num += 1
                vadv!(plasma, electricfield,
                      advection_number = vel_adv_num,
                      filtering = velocity_filtering && (vel_adv_num == 1 )) # Apply filter just once per time iteration
            end
        end

        # Save data and notify progress at checkpoints
        datasaver(plasma, time_manager, t, progress_file = progress_file)
    end
    return 0;
end
