
function (integrator::VlasovaIntegrator)(plasma::Plasma,
                                         Nt::Integer, dt::Float64,
                                         poisson!::Poisson,
                                         external_potential::Function,
                                         sadv!::SpaceAdvection,
                                         vadv!::VelocityAdvection,
                                         velocity_filtering::Bool,
                                         datasaver::DataSaver;
                                         progress_file::String)

    # Preallocated to make operations in place
    chargedensity = get_density( plasma )
    electricfield = poisson!( chargedensity )
    
    checkpoint_percent = fld(100, length(datasaver.checkpoint_axis) - 1)
    iteration_axis = (datasaver.last_iteration_saved + 1):Nt
    # TODO: Allow for merging
    notify("Entering main loop...", filename = progress_file, mode = "w")
    start_time = Dates.now()
    for t in iteration_axis     # Time loop
        time = (t-2) * dt
        pos_adv_num = 0
        vel_adv_num = 0
        for a in integrator.sequence
            if a == 'A'
                pos_adv_num += 1
                sadv!(plasma, advection_number = pos_adv_num)
                time += sadv!.coefficients[ pos_adv_num ]
                get_density!(chargedensity, plasma)
                poisson!(electricfield, chargedensity, external_potential = external_potential( time, plasma.box ) )
           else
                vel_adv_num += 1
                vadv!(plasma, electricfield,
                      advection_number = vel_adv_num,
                      filtering = velocity_filtering && (vel_adv_num == 1 )) # Apply filter just once per iteration
            end
        end
        datasaver(plasma, t)

        savedf = any( isapprox.(time, datasaver.save_distribution_times, atol = 1e-10 ) )
        savedf ? save_distribution(datasaver, plasma) : nothing
            
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
