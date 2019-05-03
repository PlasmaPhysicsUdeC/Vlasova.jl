"""
    This struct carries all te information needed to solve the Poisson equation efficiently.
        Usage:
            julia> poisson! = Poisson(plasma::Plasma)
            julia> poisson!(electricfield, chargedensity)
    where chargedensity is the charge density of the plasma and electricfield is preallocated 
"""
mutable struct Poisson          # Todo: mutable, isnt it?
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
                                   plasma.box.space_dims, flags = FFTW_flags ) # TODO: FFTW_flags
            
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
            for d in 2:plasma.box.dim_axis[end]
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

"""
    This function solves the Poisson equation to find the electric field generated
    by a charge distribution:
    
       \\nabla \\cdot \\vec E = -4\\pi \\rho
    where \\vec E is the electric field and \\rho is the charge density.

    When called with the electricfield and the chargedensity it works in place, but
    if its called with just the chargedensity, it returns the electricfield.
"""
function (poisson::Poisson)(electricfield, chargedensity )
    LinearAlgebra.mul!( poisson.fourier_density, poisson.plan, chargedensity)
    
    for d in length(electricfield)
        LinearAlgebra.ldiv!( electricfield[d], poisson.plan, poisson.integrate[d] .* poisson.fourier_density )
    end
    return 0;
end

"""
    Takes only the chargedensity and returns the electric field
"""
function (poisson::Poisson)(chargedensity)
    Nx = size( chargedensity )
    electricfield = Array{Array{Float64}}(undef, length(Nx))
    for d in 1:length(Nx)
        electricfield[d] = similar(chargedensity)
    end
    poisson( electricfield, chargedensity) # Call in-place version
    return electricfield;
end
