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
        Nx2p1, fourier_space_indices = get_rfft_dims( plasma.box.x )
        k = rfft_wavevector( plasma.box.x )

        # advection coefficients from integrator
        pos_ind = findall([ i == 'A' for i in integrator.sequence ])
        pos_coefficients = integrator.coefficients[ pos_ind ] * dt
        
        # specie coefficients in fourier space
        specie_coefficients = [ sqrt( plasma.species[s].temperature /
                                      plasma.species[s].mass ) for s in plasma.specie_axis]
        
        ktimesv = zeros(Float64, size(transformed_DF) )
        for d in plasma.box.dim_axis, i in fourier_space_indices, j in CartesianIndices(plasma.box.Nv)
            ktimesv[i, j] += k[d][i[d]] * plasma.box.v[d][j[d]]
        end
        # shift: Array (advections) of array (species) of space propagators
        shift = [ [ exp.( - 1im * pos_coeff * sp_coeff * ktimesv ) for sp_coeff in specie_coefficients ]
                  for pos_coeff in pos_coefficients ]
        # Make struct
        new( pos_coefficients,
             plan,
             transformed_DF,
             shift )
    end
end
