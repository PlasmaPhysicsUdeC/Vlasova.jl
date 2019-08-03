struct SpaceAdvection
    coefficients::Array{Float64}
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    shift::Array{Array{Array{Complex{Float64}}}}

    # Constructor from a plasma, dt and [optionally] FFT_flags
    SpaceAdvection(plasma::Plasma, integrator::VlasovaIntegrator, dt::Float64; FFTW_flags = FFTW.ESTIMATE ) = begin
        plan = FFTW.plan_rfft( copy(plasma.species[1].distribution),
                               plasma.box.space_dims, flags = FFTW_flags )
        transformed_DF = plan * plasma.species[1].distribution
        # Prepare space shift
        N2p1 = Int32.(size(transformed_DF))
        Nx2p1 = Tuple( N2p1[i] for i in 1:plasma.box.number_of_dims )
        fourier_space_axis = CartesianIndices( Nx2p1 )
        specie_coefficients = [ sqrt( plasma.species[s].temperature /
                                      plasma.species[s].mass ) for s in plasma.specie_axis]
        k = Array{Array{Float64}}(undef, plasma.box.number_of_dims)
        k[1] = rfft_wavevector( plasma.box.x[1] )
        for d in 2:plasma.box.number_of_dims
            k[d] = wavevector( plasma.box.x[d] )
        end
        pos_ind = findall([ i == 'A' for i in integrator.sequence ])
        pos_coefficients = integrator.coefficients[pos_ind] * dt
        number_of_space_advections = length( pos_coefficients )
        
        shift = Array{Array{Array{Complex{Float64}}}}(undef, number_of_space_advections)
        for a in 1:number_of_space_advections
            tmp = zeros(Float64, N2p1 )
            for d in plasma.box.dim_axis, i in fourier_space_axis, j in plasma.box.velocity_axis
                tmp[i, j] += k[d][i[d]] * plasma.box.v[d][j[d]]
            end
            shift[a] = Array{Array{Complex{Float64}}}(undef, plasma.number_of_species)
            for s in plasma.specie_axis
                shift[a][s] = exp.( - 1im * pos_coefficients[a] * specie_coefficients[s] * tmp )
            end
        end
        # Make struct
        new( pos_coefficients,
             plan,
             transformed_DF,
             shift )
    end
end
