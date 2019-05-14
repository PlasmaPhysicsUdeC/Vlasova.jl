"""
    This struct carries all te information needed to solve the Poisson equation efficiently.
        Usage:
            julia> poisson! = Poisson(plasma::Plasma)
            julia> poisson!(electricfield, chargedensity)
    where chargedensity is the charge density of the plasma and electricfield is preallocated 
"""
struct Poisson          # Todo: mutable, isnt it?
    fourier_density::Array{Complex{Float64}}
    integrate::Array{Array{Complex{Float64}}}
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
            
            integrate = Array{Array{Complex{Float64}}}(undef, plasma.box.number_of_dims)
            for d in 1:plasma.box.number_of_dims
                integrate[d] = -1im ./ k2
                for i in fourier_axis
                    integrate[d][ i ] *= k[d][ i[d] ]
                end
            end
            
            # Make struct
            new( fourier_density,
                 integrate,
                 plan )
        end
end
