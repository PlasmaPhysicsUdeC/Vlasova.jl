"""
    This struct carries all te information needed to solve the Poisson equation efficiently.
        Usage:
            julia> poisson! = Poisson(plasma::Plasma)
            julia> poisson!(electricfield, chargedensity)
    where chargedensity is the charge density of the plasma and electricfield is preallocated 
"""
struct Poisson          # Todo: mutable, isnt it?
    fourier_density::Array{Complex{Float64}}
    k::Array{Array{Float64}}
    pot2dens::Array{Complex{Float64}}
    dens2field::Array{Array{Complex{Float64}}}
    plan::FFTW.FFTWPlan

    # Construct a Poisson struct from a Plasma and [optionally] FFTW flags
    Poisson(plasma::Plasma; FFTW_flags = FFTW.ESTIMATE) =
        begin
            plan = FFTW.plan_rfft( Array{Float64}(undef, plasma.box.Nx),
                                   plasma.box.space_dims, flags = FFTW_flags )
            
            fourier_density = zeros(Complex{Float64}, Nx2p1 )

            k = rfft_wavevector( plasma.box.x )
            Nx2p1, fourier_axis = get_rfft_dims( plasma.box.x )

            # Wavevector squared
            k2 = get_k2( plasma.box )
            minus_im_over_k2 = -1im ./ k2; minus_im_over_k2[1] = 0.0

            # Wavevector repeated to match density size
            k_multidim = [[ k[d][i[d]] for i in CartesianIndices( Nx2p1 ) ]
                          for d in plasma.box.dim_axis ]

            # Transform fourier density to fourier field
            dens2field = [ k_multidim .* minus_im_over_k2
                           for d in plasma.box.dim_axis]
            
            # Make struct
            new( fourier_density,
                 k_multidim,
                 k2,
                 dens2field,
                 plan)
        end
end
