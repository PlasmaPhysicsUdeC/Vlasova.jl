mutable struct VelocityAdvection
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    wavevector::Array{Array{Float64}}
    filter::Array{Array{Float64}}
    N2p1::NTuple{N, Int32} where N # Int32 to be passed to Fortran!
    advection_coefficients::Array{Complex{Float64}, 1}
    specie_coefficients::Array{Float64, 1}

    # Constructor from a plasma, dt and [optionally] FFTW_flags
    VelocityAdvection(plasma::Plasma, dt; FFTW_flags = FFTW.ESTIMATE ) =
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
            space_coefficients, advection_coefficients = integration_coefficients( dt )
            
            specie_coefficients = [ plasma.species[s].charge / sqrt(
                plasma.species[s].temperature * plasma.species[s].mass ) for s in plasma.specie_axis ]
            
            # Make struct
            new(plan,
                transformed_DF,
                v_wavevector,
                filter,
                N2p1,
                advection_coefficients,
                specie_coefficients )
        end
end

"""
        Filter to get rid of numerical problems by artificially damping non-physical frequencies
        This function requires the vector of frequencies and returns the shape of the filter 
        in the frequency space.
    """
function anisotropic_filter(u::Array{Float64})
    u_max = maximum(u)
    return @. exp( -36*( u / u_max )^36 )
end

function (obj::VelocityAdvection)(plasma::Plasma, electricfield;
    advect_twice = false, filtering = true, advection_number = 1 )
    
    for s in 1:plasma.number_of_species
        
        # Calculate the coefficient for this specie and advection
        coefficient = - 1im * (obj.advection_coefficients[advection_number] *
                               obj.specie_coefficients[s] )
        
        # Apply advection twice?
        advect_twice ? (coefficient *= 2) : nothing
        
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



