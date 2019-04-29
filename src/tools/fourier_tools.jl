export wavevector, rfft_wavevector, rfft_derivative

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

function anisotropic_filter( array )
    return @. exp( -36 * ( array / maximum(array)) )
end
