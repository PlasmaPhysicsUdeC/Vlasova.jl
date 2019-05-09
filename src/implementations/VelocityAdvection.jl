"""
    Filter to get rid of numerical problems by artificially damping very high (non-physical) frequencies
    This function requires the vector of frequencies and returns the shape of the filter 
    in the frequency space.
"""
function anisotropic_filter(u::Array{Float64})
    u_max = maximum(u)
    return @. exp( -36*( u / u_max )^36 )
end

function (obj::VelocityAdvection)(plasma::Plasma, electricfield;
                                  filtering = false, advection_number = 1 )
    
    for s in 1:plasma.number_of_species
        
        # Calculate the coefficient for this specie and advection
        coefficient = - 1im * (obj.advection_coefficients[advection_number] *
                               obj.specie_coefficients[s] )
        
        # DF to velocity Fourier-space
        LinearAlgebra.mul!(obj.transformed_DF, obj.plan, plasma.species[s].distribution )
        
        # Apply high-freq filter
        filtering ? lowpass_velocityfilter!( obj.transformed_DF, obj.filter, obj.N2p1... ) : nothing
        
        # Apply advection
        _velocity_advection!( obj.transformed_DF, electricfield, obj.wavevector,
                              coefficient, obj.N2p1...)

        # return to real space
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
