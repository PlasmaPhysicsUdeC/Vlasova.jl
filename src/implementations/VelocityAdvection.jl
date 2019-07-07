"""
    Filter to get rid of numerical problems by artificially damping very high (non-physical) frequencies
    This function requires the vector of frequencies and returns the shape of the filter 
    in the frequency space.
"""
function anisotropic_filter(u::Array{Float64})
    u_max = maximum(u)
    return @. exp( -36*( u / u_max )^36 )
end

function (obj::VelocityAdvection)(plasma::Plasma, electricfield, grad, ecorrected;
                                  advection_number::Int64,
                                  gradient_number::Int64,
                                  is_gradient_advection::Bool,
                                  filtering::Bool)
    
    for s in 1:plasma.number_of_species
        
        # DF to velocity Fourier-space
        LinearAlgebra.mul!(obj.transformed_DF, obj.plan, plasma.species[s].distribution )
        
        # Apply high-freq filter
        filtering ? lowpass_velocityfilter!( obj.transformed_DF, obj.filter, obj.N2p1... ) : nothing

        # Calculate the coefficient for this specie and advection
        coefficient = - 1im * (obj.advection_coefficients[advection_number] *
                               obj.specie_coefficients[s] )

        # Apply advection
        if is_gradient_advection
            for i in 1:plasma.box.number_of_dims # TODO: Get the k/k2 to poisson to obtain grad in d>1
                @. ecorrected[i] = electricfield[i] + obj.specie_coefficients[s] * obj.gradient_coefficients[ gradient_number ] * grad[i]
            end
            _velocity_advection!( obj.transformed_DF, ecorrected, obj.wavevector,
                                  coefficient, obj.N2p1...)
        else
            _velocity_advection!( obj.transformed_DF, electricfield, obj.wavevector,
                                  coefficient, obj.N2p1...)
        end 

        # Return to real space
        LinearAlgebra.ldiv!( plasma.species[s].distribution, obj.plan, obj.transformed_DF )
        
    end
    return 0;
end


# 1D
function lowpass_velocityfilter!(transformed_DF, filter, Nx::Int32, Nvx::Int32 )
    ccall((:velocity_filter1d_, "velocity_filter.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Array{Float64, 1}}, Ptr{Array{Complex{Float64}, 2}}),
          Ref(Nx), Ref(Nvx), filter[1], transformed_DF)
    return 0;
end

# 2D
function lowpass_velocityfilter!(transformed_DF, filter, Nx::Int32, Ny::Int32, Nvx::Int32, Nvy::Int32)
    ccall((:velocity_filter2d_, "velocity_filter.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Array{Float64, 1}}, Ptr{Array{Float64, 1}}, Ptr{Array{Complex{Float64}, 4}}),
          Ref(Nx), Ref(Ny), Ref(Nvx), Ref(Nvy), filter[1], filter[2], transformed_DF)
    return 0;
end

# 1D
function _velocity_advection!(transformed_DF, electricfield, wavevector, coef, Nx::Int32, Nvx::Int32)
    ccall((:velocity_advection1d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Complex{Float64}}, Ptr{Array{Float64, 1}}, Ptr{Array{Float64, 1}}, Ptr{Array{Complex{Float64}, 2}}),
          Ref(Nx), Ref(Nvx), Ref(coef), electricfield[1], wavevector[1], transformed_DF)
    return 0;
end

# 2D
function _velocity_advection!(transformed_DF, electricfield, wavevector, coef, Nx::Int32, Ny::Int32, Nvx::Int32, Nvy::Int32)
    ccall((:velocity_advection2d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex{Float64}}, Ptr{Array{Float64, 2}},
           Ptr{Array{Float64, 2}}, Ptr{Array{Float64, 1}}, Ptr{Array{Float64, 1}}, Ptr{Array{Complex{Float64}, 4}}),
          Ref(Nx), Ref(Ny), Ref(Nvx), Ref(Nvy), Ref(coef), electricfield[1], electricfield[2],
          wavevector[1], wavevector[2], transformed_DF)
    return 0;
end
