mutable struct SpaceAdvection
    plan::FFTW.FFTWPlan
    transformed_DF::Array{Complex{Float64}}
    shift::Array{Array{Array{Complex{Float64}}}}
    N2p1::NTuple{N, Int32} where N # Int32 to be passed to Fortran!

    # Constructor from a plasma, dt and [optionally] FFT_flags
    SpaceAdvection(plasma::Plasma, dt::Float64; FFTW_flags = FFTW.ESTIMATE ) = begin
        plan = FFTW.plan_rfft( copy(plasma.species[1].distribution),
                               plasma.box.space_dims, flags = FFTW_flags )
        transformed_DF = plan * plasma.species[1].distribution
        # Prepare space shift
        N2p1 = Int32.(size(transformed_DF))
        Nx2p1 = Tuple( N2p1[i] for i in 1:plasma.box.number_of_dims )
        fourier_space_axis = CartesianIndices( Nx2p1 )
        specie_coefficients = [ sqrt( plasma.species[s].temperature *
                                      plasma.species[s].mass ) for s in plasma.specie_axis]
        k = Array{Array{Float64}}(undef, plasma.box.number_of_dims)
        k[1] = rfft_wavevector( plasma.box.x[1] )
        for d in 2:plasma.box.number_of_dims
            k[d] = wavevector( plasma.box.x[d] )
        end
        pos_coefficients, vel_coefficients = integration_coefficients( dt )
        number_of_space_advections = length( pos_coefficients )
        
        shift = Array{Array{Array{Complex{Float64}}}}(undef, number_of_space_advections)
        for a in 1:number_of_space_advections
            tmp = zeros(Float64, N2p1 )
            for d in plasma.box.dim_axis, i in fourier_space_axis, j in plasma.box.velocity_axis
                tmp[i, j] += k[d][i[d]] * plasma.box.v[d][j[d]]
            end
            shift[a] = Array{Array{Complex{Float64}}}(undef, plasma.number_of_species)
            for s in plasma.specie_axis
                shift[a][s] = exp.( tmp * pos_coefficients[a] * specie_coefficients[s] )
            end
        end
        # Make struct
        new( plan,
             transformed_DF,
             shift,
             N2p1 )
    end
end

"""
    Advects a plasma in space efficiently using the information stored in the
    object space_advection of type SpaceAdvection.
"""
function (space_advection::SpaceAdvection)(plasma::Plasma; advection_number::Integer = 1)
    for s in plasma.specie_axis
        LinearAlgebra.mul!(space_advection.transformed_DF,
                           space_advection.plan, plasma.species[s].distribution)
        # Here happens the magic!
        _space_advection!( space_advection.transformed_DF,
                           space_advection.shift[advection_number][s],
                           space_advection.N2p1...)
        
        LinearAlgebra.ldiv!(plasma.species[s].distribution,
                            space_advection.plan, space_advection.transformed_DF)
    end
    return 0;
end

# 1D - Fortran
function _space_advection!(transformed_DF, shift, Nx::Int32, Nvx::Int32)
    ccall((:space_advection1d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Array{Complex{Float64}, 2}}, Ptr{Array{Complex{Float64}, 2}}),
          Ref(Nx), Ref(Nvx), shift, transformed_DF)
end

# 2D - Fortran
function _space_advection!(transformed_DF, shift, Nx::Int32, Ny::Int32, Nvx::Int32, Nvy::Int32)
    ccall((:space_advection2d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Array{Complex{Float64}, 4}}, Ptr{Array{Complex{Float64}, 4}}),
          Ref(Nx), Ref(Ny), Ref(Nvx), Ref(Nvy), shift, transformed_DF)
    return 0;
end
