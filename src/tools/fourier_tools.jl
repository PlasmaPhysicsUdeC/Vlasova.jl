export wavevector,
    rfft_wavevector,
    get_rfft_dims

"""
```julia
get_k2( box::Box )
```

Obtain the squared space-wavevector,

``k^2 = k_1^2 + k_2^2 + ... k_n^2``

where ``k_1`` is an rfft-wavevector (non-negative half of the whole wavevector).
"""
function get_k2( box::Box )
    Nx2p1, fourier_axis = get_rfft_dims( box.x )

    k = rfft_wavevector( box.x )

    k2 = zeros( Nx2p1 )
    for d in 1:box.number_of_dims, i in fourier_axis
        k2[i] += ( k[d][ i[d] ] )^2
    end

    return k2
end

"""
```julia
anisotropic_filter(box::Box)
```

Filter to get rid of numerical problems by artificially damping very high (non-physical) frequencies.

# Notes
* This function is not performance-critical since it is executed just once when [`integrate!`](@ref) is prepared.
* The result of this function is multiplied element-wise with the rfft-velocity-transformed distribution
function(s) to damp the high frequencies.
* Please take into consideration that using this kind of filtering is equivalent to adding hyperviscosity
  to the Vlasov equation.
* If you want to disable filtering, just set the Vlasova variable [`velocity_filtering`](@ref)
  to `false` using the [`@vlasova`](@ref) macro, but take into account that recurrence may appear.
* If you know what you are doing, you can redefine this function to filter in an alternate way.
"""
function anisotropic_filter(box::Box)
    u = rfft_wavevector( box.v )
    Nv2p1, fourier_indices = get_rfft_dims( box.v )

    filter = ones( Nv2p1 )
    for d in 1:box.number_of_dims
        umax = maximum(u[d])
        filter1d = @. exp( -36*( u[d] / umax )^36 )
        for i in fourier_indices
            filter[i] *= filter1d[ i[d] ]
        end
    end

    return filter
end

# ===== Exported functions start here =====

"""
```julia
wavevector(vector::Array{Float64, 1})
```

Obtain the Fourier-conjugate of a 1-dimensional array.

# Examples

```julia
julia> time_axis = collect( 0:dt:final_time );

julia> frequency = wavevector( time_axis );

```
"""
function wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N;
    wavevector = Array(1:N) .- (N/2 +1);

    return FFTW.fftshift( wavevector ) * 2pi/length;
end


"""
```julia
rfft_wavevector(vector::Array{Float64, 1})
```

Obtain the rfft-Fourier-conjugate of a 1-dimensional array.

# Examples

```julia
julia> kx = rfft_wavevector( box.x[1] );
```
"""
function rfft_wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N

    return Array( 0:div(N,2) ) * 2pi/length;
end

"""
```julia
rfft_wavevector(vector::Array{Array{Float64, 1}, 1} )
```

Obtain the Fourier-conjugate of the space or velocity variables
using the usual conventions of Vlasova.

# Examples

```julia
julia> k = rfft_wavevector( box.x );

```
"""
function rfft_wavevector(vector::Array{Array{Float64, 1}, 1} )
    number_of_dims = length(vector)

    k = Array{Array{Float64, 1}}( undef, number_of_dims )
    k[1] = rfft_wavevector( vector[1] )
    for d in 2:number_of_dims
        k[d] = wavevector( vector[d] )
    end

    return k
end

"""
```julia
get_rfft_dims(A::Array{Float64}; transformed_dims)
```

Obtain the `size` and `CartesianIndices` of an array after it is rfft-transformed.
"""
function get_rfft_dims(A::Array{Float64}; transformed_dims)
    N = size( A )
    @assert maximum( transformed_dims ) < length( N) "all transformed_dims must exist in the array"

    rdim = transformed_dims[1]
    N2p1 = Tuple( (i == rdim ) ? div( N[rdim], 2 ) + 1 : N[i] for i in 1:size(N, 1) )
    fourier_axis = CartesianIndices( N2p1 )

    return N2p1, fourier_axis
end

"""
```julia
get_rfft_dims(x::Array{Array{Float64, 1}})
```
Obtain the `size` and `CartesianIndices` of the Fourier-conjugate variable of the input.

# Important
* This extension of the function was developed to be used only over `box.x` and `box.v`
  inside the inner functions of Vlasova. Any other use is not encouraged.
"""
function get_rfft_dims(x::Array{Array{Float64, 1}})
    sz = Tuple( size(a, 1) for a in x )
    Nx2p1 = Tuple( (i == 1) ? div( sz[i], 2 ) + 1 : sz[i] for i in 1:size(sz, 1) )

    fourier_axis = CartesianIndices( Nx2p1 )

    return Nx2p1, fourier_axis
end


"""
```julia
get_rfft_dims(box::Box)
```
Obtain the `size` and `CartesianIndices` of the whole Fourier space of a [`Box`](@ref),
including the space and velocity dimensions.
"""
function get_rfft_dims(box::Box)

    N2p1 = Tuple( (i == 1) ? div( box.N[1], 2 ) + 1 : box.N[i] for i in 1:box.number_of_dims )
    fourier_axis = CartesianIndices( N2p1 )

    return N2p1, fourier_axis
end
