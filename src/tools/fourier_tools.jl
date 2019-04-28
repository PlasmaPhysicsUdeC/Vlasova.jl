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

function rfft_derivative(x::Array{Float64}, k::Array{Complex{Float64}}; ord::Integer = 1, dims = 1) # 2DMax __Under DEV__
    # Derivative of real arrays by use of rfft
    # This function also performs integrals (ord < 0)
    @assert (dims == 1) || (dims == 2) "rfftDerivative is only defined along dims 1 or 2"
    kk = (1im*k).^ord; kk[1] = 0.0im;
    if dims == 2
        kk = kk'
    end
    return FFTW.irfft(rfft(x, [dims]).*kk , size(x, dims), [dims] )
end
