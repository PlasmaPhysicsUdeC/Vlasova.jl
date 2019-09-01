"""
```julia
Box(;Nx, Nv, Lx, vmin, vmax)
```

Create an object to store all the information about the physical space under which a [`Vlasova.Plasma`](@ref) lives.

# Notes
* All arguments to create a box are keywords.
* `Nx` and `Nv` must be an `Integer`, a `Tuple{Integer}` or an `Array{Integer}`.
* `Lx`, `vmin` and `vmax` must be a `Real`, a `Tuple{Real}` or an `Array{Real}`.
* All of the keywords given must have the same number of elements, or an error will be thrown.
* All `Real`s will be converted to `Float64` and all `Integer`s will be converted to `Int64`.
"""
struct Box
    # Fundamental quantities
    Nx::NTuple{N, Int64} where N
    Nv::NTuple{N, Int64} where N
    Lx::NTuple{N, Float64} where N
    vmin::NTuple{N, Float64} where N
    vmax::NTuple{N, Float64} where N
    # Derived quantities
    ## Dims
    N::NTuple{N, Int64} where N
    number_of_dims::Int64
    space_dims::NTuple{N, Int64} where N
    velocity_dims::NTuple{N, Int64} where N
    ## Axes
    dim_axis::Base.OneTo{Int64}
    space_axes::NTuple{N, Base.OneTo{Int64}} where N
    velocity_axes::NTuple{N, Base.OneTo{Int64}} where N
    distribution_axes::NTuple{N, Base.OneTo{Int64}} where N
    ## Indices
    space_indices::CartesianIndices{N, NTuple{N, Base.OneTo{Int64}}} where N
    velocity_indices::CartesianIndices{N, NTuple{N, Base.OneTo{Int64}}} where N
    ## Position and velocity
    dx::NTuple{N, Float64} where N
    dv::NTuple{N, Float64} where N
    x::Array{Array{Float64, 1}}
    v::Array{Array{Float64, 1}}

    # Constructors

    # Construct with the fundamental parameters only
    _Box(Nx::NTuple{N, Int64} where N,
         Nv::NTuple{N, Int64} where N,
         Lx::NTuple{N, Float64} where N,
         vmin::NTuple{N, Float64} where N,
         vmax::NTuple{N, Float64} where N ) = begin
             # Do sizes match?
             sizes_match = ( length( Nx ) == length( Nv ) == length( Lx )
                             == length( vmin ) == length( vmax ) )
             @assert sizes_match "Sizes of Nx, Nv, Lx, vmin and vmax do not match"
             # Construct derived quantities
             ## Dims
             N = Tuple((Nx..., Nv...))
             number_of_dims = length( Nx )
             space_dims = Tuple(1:number_of_dims)
             velocity_dims = space_dims .+ number_of_dims
             ## Axes
             dim_axis = Base.OneTo(number_of_dims)
             space_axes = Base.OneTo.( Nx )
             velocity_axes = Base.OneTo.( Nv )
             distribution_axes = Base.OneTo.( N )
             ## Indices
             space_indices = CartesianIndices( space_axes )
             velocity_indices = CartesianIndices( velocity_axes )
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
             new(Nx,
                 Nv,
                 Lx,
                 vmin,
                 vmax,
                 N,
                 number_of_dims,
                 space_dims,
                 velocity_dims,
                 dim_axis,
                 space_axes,
                 velocity_axes,
                 distribution_axes,
                 space_indices,
                 velocity_indices,
                 dx,
                 dv,
                 x,
                 v )

         end

    # Convert (elements, Arrays or Tuples of) Integers/Reals into Tuples of Int64/Float64
    Box(args...
        ;Nx::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Integer,
        Nv::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Integer,
        Lx::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real,
        vmin::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real,
        vmax::Union{T, Array{T, 1}, NTuple{N, T} where N } where T <: Real
        ) = begin
            @assert args == () "Box should not be called using non-keyword arguments."
            _Box(Int64.(Tuple(Nx)),
                 Int64.(Tuple(Nv)),
                 Float64.(Tuple(Lx)),
                 Float64.(Tuple(vmin)),
                 Float64.(Tuple(vmax))
                 )
        end
end

# Display box in the REPL
function Base.display(b::Box)
    d = "
    ---
    $(b.number_of_dims)-dimensional Vlasova Box.
    ---

    Nx = $(b.Nx)
    Nv = $(b.Nv)
    Lx = $(b.Lx)
    vmin = $(b.vmin)
    vmax = $(b.vmax)
"
    print(d)
end
