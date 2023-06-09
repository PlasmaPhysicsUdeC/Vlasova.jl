struct VelocityAdvection
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    filter::Array{Float64}
    advection_coefficients::Array{Float64}
    gradient_coefficients::Array{Float64}

    # Constructor from a plasma, dt and [optionally] FFTW_flags
    VelocityAdvection(plasma::Plasma, integrator::VlasovaIntegrator, dt; FFTW_flags = FFTW.ESTIMATE ) =
        begin
            # Plan
            plan = FFTW.plan_rfft( copy(plasma.species[1].distribution),
                                   plasma.box.velocity_dims, flags = FFTW_flags )

            # transformed_DF
            transformed_DF = plan * plasma.species[1].distribution

            # Filter
            filter = anisotropic_filter( plasma.box )

            # Coefficients
            vel_ind = findall([ i in "BC" for i in integrator.sequence ])
            grad_ind = findall([ i in "C" for i in integrator.sequence ])

            # Coefficients depend upon specie and advection number

            # delta v_wavevector
            delta_p =  [ 2pi / (plasma.box.vmax[d] - plasma.box.vmin[d])
                         for d in 1:plasma.box.number_of_dims ]

            # Specie coefficients
            adv_spc_coefficients = [ plasma.species[s].charge / sqrt(
                plasma.species[s].temperature * plasma.species[s].mass )
                                     for s in plasma.specie_axis ]

            grad_spc_coefficients = [ plasma.species[s].charge^2 / sqrt(
                plasma.species[s].mass^3 * plasma.species[s].temperature )
                                      for s in plasma.specie_axis ]

            # Total coefficients
            advection_coefficients = dt * [ adv_coeff * spc_coeff * dp
                                            for dp in delta_p,
                                            spc_coeff in adv_spc_coefficients,
                                            adv_coeff in integrator.coefficients[vel_ind] ]

            gradient_coefficients = dt^3 * [ grad_coeff * spc_coeff * dp
                                             for dp in delta_p,
                                             spc_coeff in grad_spc_coefficients,
                                             grad_coeff in ( integrator.gradient_coefficients .*
                                                             integrator.coefficients[grad_ind] ) ]

            # Make struct
            new(plan,
                transformed_DF,
                filter,
                advection_coefficients,
                gradient_coefficients )
        end
end
