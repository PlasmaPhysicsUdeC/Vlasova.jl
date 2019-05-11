export wavevector, rfft_wavevector, anisotropic_filter

function get_rfft_dims(box::Box)
    Nx2p1 = Tuple( (i == 1) ? fld( box.Nx[1], 2 ) + 1 : box.Nx[i] for i in 1:box.number_of_dims )
    fourier_axis = CartesianIndices( Nx2p1 )
    
    return Nx2p1, fourier_axis
end

function wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N;
    wavevector = Array(1:N) .- (N/2 +1);
    return FFTW.fftshift( wavevector ) * 2pi/length;
end

function rfft_wavevector(vector::Array{Float64, 1})
    N = size(vector, 1)
    length = (vector[2] - vector[1] )*N
    return Array( 0:fld(N,2) ) * 2pi/length;
end

"""
    Obtain the wavevector of an array such as position or velocity
    with the form used in Vlasova.jl
"""
function rfft_wavevector(vector::Array{Array{Float64, 1}, 1} )
    k = Array{Array{Float64, 1}}( undef, length(vector) )
    k[1] = rfft_wavevector( vector[1] )
    for d in 2:box.number_of_dims
        k[d] = wavevector( vector[d] )
    end
    return k
end

function anisotropic_filter( array )
    return @. exp( -36 * ( array / maximum(array)) )
end
