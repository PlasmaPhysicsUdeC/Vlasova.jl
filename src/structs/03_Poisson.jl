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
            Nx2p1 = Tuple( i == 1 ? fld(plasma.box.Nx[i], 2)+1 : plasma.box.Nx[i]
                           for i in 1:length(plasma.box.Nx) )
            fourier_axis = CartesianIndices( Nx2p1 )
            
            plan = FFTW.plan_rfft( Array{Float64}(undef, plasma.box.Nx),
                                   plasma.box.space_dims, flags = FFTW_flags )
            
            fourier_density = zeros(Complex{Float64}, Nx2p1 )
            
            k = rfft_wavevector( plasma.box.x )
            k2 = get_k2( plasma.box )
            k2[1] = Inf  # So that the inverse yields 0.0. Ensure quasineutrality
            
            dens2field = Array{Array{Complex{Float64}}}(undef, plasma.box.number_of_dims)
            k_multidim = Array{Array{Float64}}(undef, plasma.box.number_of_dims)
            for d in 1:plasma.box.number_of_dims
                k_multidim[d] = ones( Nx2p1... )
                dens2field[d] = -1im ./ k2
                for i in fourier_axis
                    k_multidim[d][ i ] *= k[d][ i[d] ]
                    dens2field[d][ i ] *= k[d][ i[d] ]
                end
            end
            k2[1] = 0.0im
            
            # Make struct
            new( fourier_density,
                 k_multidim,
                 k2,
                 dens2field,
                 plan)
        end
end
