export wavevector,
    get_rfft_dims,
    rfft_wavevector,
    anisotropic_filter

"""
    Obtain the size of an array after it is rfft-transformed
    The function also outputs the CartesianIndices of the transformed object
"""
function get_rfft_dims(A::Array{Float64}; transformed_dims)
    N = size( A )
    @assert maximum( transformed_dims ) < length( N) "all transformed_dims must exist in the array"

    rdim = transformed_dims[1]
    N2p1 = Tuple( (i == rdim ) ? div( N[rdim], 2 ) + 1 : N[i] for i in 1:size(N, 1) )
    fourier_axis = CartesianIndices( N2p1 )
    
    return N2p1, fourier_axis
end

function get_rfft_dims(x::Array{Array{Float64, 1}})
    sz = Tuple( size(a, 1) for a in x )
    Nx2p1 = Tuple( (i == 1) ? div( sz[i], 2 ) + 1 : sz[i] for i in 1:size(sz, 1) )

    fourier_axis = CartesianIndices( Nx2p1 )
    
    return Nx2p1, fourier_axis
end

function get_rfft_dims(box::Box)
    
    N2p1 = Tuple( (i == 1) ? div( box.N[1], 2 ) + 1 : box.N[i] for i in 1:box.number_of_dims )
    fourier_axis = CartesianIndices( N2p1 )
    
    return N2p1, fourier_axis
end

"""
    Obtain the squared wavevector,

    `k^2 = k_1^2 + k_2^2 + ... k_n^2`

    where k_1 is a half of the first wavevector as it would be used on a real DFT.
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

# Exported functions start here

function wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N;
    wavevector = Array(1:N) .- (N/2 +1);
    
    return FFTW.fftshift( wavevector ) * 2pi/length;
end

function rfft_wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N
    
    return Array( 0:div(N,2) ) * 2pi/length;
end

"""
    Obtain the wavevector of an array such as position or velocity
    with the form used in Vlasova.jl
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
        Filter to get rid of numerical problems by artificially damping very high (non-physical) frequencies
        This function requires the vector of frequencies and returns the shape of the filter 
        in the frequency space.
"""
function anisotropic_filter(box::Box)
    wavevector = rfft_wavevector( box.v )
    Nv2p1, fourier_indices = get_rfft_dims( box.v )

    filter = ones( Nv2p1 )
    for d in box.number_of_dims
        filter1d = _anisotropic_filter( wavevector[d] )
        for i in fourier_indices
            filter[i] *= filter1d[ i[d] ]
        end
    end

    return filter
end

"""
    Define the shape of the anisotropic filter along each dimension.

    It is unexported and any extension to the anisotropic filter is recommended to be performed on
    the function anisotropic_filter(::Box) rather than this one, since that is the one called when 
    a VelocityAdvection is built.
"""
function _anisotropic_filter(u::Array{Float64}) # 1-dimensional
    u_max = maximum(u)

    return @. exp( -36*( u / u_max )^36 )
end
