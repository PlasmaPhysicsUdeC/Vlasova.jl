struct VelocityAdvection
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    wavevector::Array{Array{Float64}}
    filter::Array{Array{Float64}}
    N2p1::NTuple{N, Int32} where N # Int32 to be passed to Fortran!
    advection_coefficients::Array{Complex{Float64}, 1}
    gradient_coefficients::Array{Complex{Float64}, 1}
    specie_coefficients::Array{Float64, 1}

    # Constructor from a plasma, dt and [optionally] FFTW_flags
    VelocityAdvection(plasma::Plasma, integrator::VlasovaIntegrator, dt; FFTW_flags = FFTW.ESTIMATE ) =
        begin
            plan = FFTW.plan_rfft( copy(plasma.species[1].distribution),
                                   plasma.box.velocity_dims, flags = FFTW_flags )
            
            transformed_DF = plan * plasma.species[1].distribution
            N2p1 = Int32.( size(transformed_DF) )
            
            v_wavevector = Array{Array{Float64}}(undef, plasma.box.number_of_dims)

            v_wavevector[1] = rfft_wavevector( plasma.box.v[1] )            
            for d in 2:plasma.box.number_of_dims
                v_wavevector[d] = wavevector( plasma.box.v[d] )
            end
            
            # Filter
            filter = Array{Array{Float64}}(undef, plasma.box.number_of_dims)
            for d in 1:plasma.box.number_of_dims
                filter[d] = anisotropic_filter( v_wavevector[d] )
            end

            # Coefficients
            ## advection
            vel_ind = findall([ i in "BC" for i in integrator.sequence ])
            advection_coefficients = integrator.coefficients[vel_ind] * dt
            gradient_coefficients = integrator.gradient_coefficients * dt^2
            ## specie
            specie_coefficients = [ plasma.species[s].charge / sqrt(
                plasma.species[s].temperature * plasma.species[s].mass ) for s in plasma.specie_axis ]
            
            # Make struct
            new(plan,
                transformed_DF,
                v_wavevector,
                filter,
                N2p1,
                advection_coefficients,
                gradient_coefficients,
                specie_coefficients )
        end
end
