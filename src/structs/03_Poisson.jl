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
            
            k = Array{Array{Float64, 1}}(undef, plasma.box.number_of_dims)

            k[1] = rfft_wavevector( plasma.box.x[1] )
            for d in 2:plasma.box.number_of_dims
                k[d] = wavevector( plasma.box.x[d] )
            end
            
            k2 = zeros( Nx2p1 )
            for d in 1:plasma.box.number_of_dims, i in fourier_axis
                k2[i] += ( k[d][ i[d] ] )^2
            end
            k2[1] = Inf  # So that the inverse yields 0.0. Ensure quasineutrality
            
            integrate = Array{Array{Complex{Float64}}}(undef, plasma.box.number_of_dims)
            integrate[1] = -1im ./ k2
            for d in 2:plasma.box.number_of_dims
                integrate[d] = integrate[1]
            end
            
            for d in plasma.box.dim_axis, i in fourier_axis
                integrate[d][ i ] *= k[d][ i[d] ]
            end
            # Make struct
            new( fourier_density,
                 integrate,
                 plan )
        end
end
