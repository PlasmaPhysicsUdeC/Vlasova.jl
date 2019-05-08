# To avoid closures, parameters will be passed between functions in a container

"""
    Store the information about the physical Box inside which the simulation runs. #TODO
"""
struct Box
    # Fundamental quantities
    simulation_name::String
    Nx::NTuple{N, Int64} where N
    Nv::NTuple{N, Int64} where N
    Lx::NTuple{N, Float64} where N
    vmin::NTuple{N, Float64} where N
    vmax::NTuple{N, Float64} where N
    # Derived quntities
    ## Dims
    N::NTuple{N, Int64} where N
    number_of_dims::Int64
    space_dims::NTuple{N, Int64} where N
    velocity_dims::NTuple{N, Int64} where N
    ## Axes
    dim_axis::Base.OneTo{Int64}
    space_axis::CartesianIndices{M, NTuple{M, Base.OneTo{Int64}}} where M
    velocity_axis::CartesianIndices{M, NTuple{M, Base.OneTo{Int64}}} where M
    distribution_axis::CartesianIndices{M, NTuple{M, Base.OneTo{Int64}}} where M
    ## Position and velocity
    dx::NTuple{N, Float64} where N
    dv::NTuple{N, Float64} where N
    x::Array{Array{Float64, 1}}
    v::Array{Array{Float64, 1}}

    # Constructors

    # Construct with the fundamental parameters only
    Box(name::String,
        Nx::NTuple{N, Int64} where N,
        Nv::NTuple{N, Int64} where N,
        Lx::NTuple{N, Float64} where N,
        vmin::NTuple{N, Float64} where N,
        vmax::NTuple{N, Float64} where N ) = begin
            # Do sizes match?
            if ( length( Nx ) == length( Nv ) == length( Lx )
                 == length( vmin ) == length( vmax ) )
                # Construct derived quantities
                ## Dims
                N = Tuple((Nx..., Nv...))
                number_of_dims = length( Nx )
                space_dims = Tuple(1:number_of_dims)
                velocity_dims = space_dims .+ number_of_dims
                ## Axes
                dim_axis = Base.OneTo(number_of_dims)
                space_axis = CartesianIndices( Nx )
                velocity_axis = CartesianIndices( Nv )
                distribution_axis = CartesianIndices( Tuple((Nx..., Nv...)) )
                ## Position and velocity
                Lv = @. vmax .- vmin
                dx = @. Lx / Nx
                dv = @. Lv / Nv
                x = Array{Array{Float64, 1}}(undef, number_of_dims)
                v = Array{Array{Float64, 1}}(undef, number_of_dims)
                for d in dim_axis
                    x[d] = [ i * dx[d] for i in 0:(Nx[d]-1) ]
                    v[d] = [ vmin[d] + i*dv[d] for i in 0:(Nv[d]-1) ]
                end 
                # Make the struct 
                new(name, 
                    Nx,
                    Nv,
                    Lx,
                    vmin,
                    vmax,
                    N,
                    number_of_dims,
                    space_dims,
                    velocity_dims,
                    dim_axis,
                    space_axis,
                    velocity_axis,
                    distribution_axis,
                    dx,
                    dv,
                    x,
                    v )
            else
                error("Sizes of Nx, Nv, Lx, vmin and vmax do not match")
            end
        end
    
    # Convert (elements, Arrays or Tuples of) Integers/Reals into Tuples of Int64/Float64 (1D case)
    Box(name::String,
        Nx::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Integer,
        Nv::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Integer,
        Lx::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real,
        vmin::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real,
        vmax::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real
        ) = Box(name,
                Int64.(Tuple(Nx)),
                Int64.(Tuple(Nv)),
                Float64.(Tuple(Lx)),
                Float64.(Tuple(vmin)),
                Float64.(Tuple(vmax))
                )
end
