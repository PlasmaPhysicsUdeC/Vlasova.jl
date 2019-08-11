struct VelocityAdvection
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    wavevector::Array{Array{Float64}}
    filter::Array{Float64}
    advection_coefficients::Array{Array{Float64, 1}, 1}
    gradient_coefficients::Array{Array{Float64, 1}, 1}

    # Constructor from a plasma, dt and [optionally] FFTW_flags
    VelocityAdvection(plasma::Plasma, integrator::VlasovaIntegrator, dt; FFTW_flags = FFTW.ESTIMATE ) =
        begin
            plan = FFTW.plan_rfft( copy(plasma.species[1].distribution),
                                   plasma.box.velocity_dims, flags = FFTW_flags )
            
            transformed_DF = plan * plasma.species[1].distribution
            
            wavevector = rfft_wavevector( plasma.box.v )
            Nv2p1, fourier_indices = get_rfft_dims( plasma.box.v )
            
            # Filter
            filter = anisotropic_filter( plasma.box )

            # Coefficients
            vel_ind = findall([ i in "BC" for i in integrator.sequence ])

            adv_coeffs = integrator.coefficients[vel_ind] * dt

            # TODO: This is not working but idk why. Its better when amplified by 1/dt^3
            grad_coeffs = ( integrator.gradient_coefficients  ) * dt^2

            # Coefficients depend upon specie and advection: Array{Array{Float64}
            # Specie coefficients
            adv_spc_coefficients = [ plasma.species[s].charge / sqrt(
                plasma.species[s].temperature * plasma.species[s].mass )
                                     for s in plasma.specie_axis ]
            
            grad_spc_coefficients = [ plasma.species[s].charge / plasma.species[s].mass
                                      for s in plasma.specie_axis ]
            
            # Total coefficients
            advection_coefficients = [[ adv_coeff * spc_coeff
                                        for spc_coeff in adv_spc_coefficients ]
                                      for adv_coeff in adv_coeffs ]
            
            gradient_coefficients = [[ grad_coeff * spc_coeff
                                       for spc_coeff in grad_spc_coefficients ]
                                     for grad_coeff in grad_coeffs ]
            
            # Make struct
            new(plan,
                transformed_DF,
                wavevector,
                filter,
                advection_coefficients,
                gradient_coefficients )
        end
end
