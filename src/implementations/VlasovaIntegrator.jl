
function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                         Nt::Integer, dt::Float64,
                                         poisson!::Poisson,
                                         sadv!::SpaceAdvection,
                                         vadv!::VelocityAdvection,
                                         velocity_filtering::Bool,
                                         datasaver::DataSaver;
                                         progress_file::String = "/")

    # Preallocated to make operations in place
    chargedensity = get_density( plasma )
    electricfield = poisson!( chargedensity )
    
    checkpoint_percent = fld(100, length(datasaver.checkpoint_axis) - 1)
    iteration_axis = (datasaver.last_iteration_saved + 1):Nt
    time_axis = (0:Nt-1)*dt
    
    # TODO: Allow for merging
    # if integrator.merge_last_advection
    #     if integrator.sequence[1] == 'A'
    #         sadv!(plasma, advection_number = 1)
    #     else
    #         vadv!(plasma,electricfield,
    #               advection_number = 1)
    #     end
    # end
    notify("Entering main loop...", filename = progress_file, mode = "w")
    start_time = Dates.now()
    for t in iteration_axis     # Time loop
        pos_adv_num = 0
        vel_adv_num = 0
        for a in integrator.sequence
            if a == 'A'
                pos_adv_num += 1
                sadv!(plasma, advection_number = pos_adv_num)
                get_density!(chargedensity, plasma)
                poisson!(electricfield, chargedensity)
           else
                vel_adv_num += 1
                vadv!(plasma, electricfield,
                      advection_number = vel_adv_num,
                      filtering = velocity_filtering && (vel_adv_num == 1 )) # Apply filter just once per iteration
            end
        end
        datasaver(plasma, t)
        if t in datasaver.checkpoint_axis
            save_to_disk(datasaver, plasma, t)

            # Inform progress
            accomplished = (datasaver.checkpoints_reached - 1) * checkpoint_percent
            elapsed = round(Dates.now() - start_time, Dates.Second )
            notify("\t$accomplished% accomplished in $elapsed", filename = progress_file) 
        end
    end
    return 0;
end
