export reducedims,
    outer,
    ⊗,
    cosine_perturbation1d,
    center,
    center!

"""
```julia
reducedims(f::Function, A::Array; dims = (0))
```

Perform a reduction function, `f`, over an array, `A`, along the dimension(s), `dims`, and drop the dimensions reduced.

In general,
```julia
reducedims( f, A, dims = dims)
```
is equivalent to
```julia
dropdims( f(A, dims = dims), dims = dims )
```
"""
function reducedims(f::Function, A::Array; dims = (0))
    if dims == (0)
        return f(A)
    else
        return dropdims(f(A, dims = dims), dims = dims)
    end
end

"""
```julia
outer(A::Array{TA}, B::Array{TB}, f::Function = *) where TA where TB
```
Generalized external operation between two arrays, `A` and `B`.

In general, `outer(A, B, f)` returns `R` such that

``R_{ij} = f(A_i, B_j) ``

# Notes
* If `f` is not specified, `outer` corresponds to the external product.
* This function is not optimized to be fast, but only to offer a comfortable syntax.
  It should not be used in performance-critical parts of the code.

# Examples
```jldoctest; setup = :(using Vlasova)
julia> A = [1.0, 2.0, 3.0];

julia> B = [0.0, 1.0];

julia> outer(A, B, *)
3×2 Array{Float64,2}:
 0.0  1.0
 0.0  2.0
 0.0  3.0

julia> outer(A, B, -)
3×2 Array{Float64,2}:
 1.0  0.0
 2.0  1.0
 3.0  2.0
```
"""
function outer(A::Array{TA},
               B::Array{TB},
               f::Function = *) where TA where TB
    # Sizes
    sA = size(A); sB = size(B)
    # Result
    R = Array{ promote_type(TA, TB) }(undef, sA..., sB...)
    @inbounds for j in CartesianIndices(sB)
        @inbounds for i in CartesianIndices(sA)
            R[i, j] = f( A[i], B[j] )
        end
    end
    return R
end


"""
```julia
A ⊗ B
```

Infix form of the external product between two arrays.
See [`Vlasova.outer`](@ref).

# Examples
```jldoctest; setup = :(using Vlasova)
julia> A = [1, 2, 3];

julia> B = [0, 1];

julia> A ⊗ B
3×2 Array{Int64,2}:
 0  1
 0  2
 0  3
```
"""
⊗ = Vlasova.outer

"""
```julia
cosine_perturbation1d(box::Box; modes, amplitudes, dim = 1)
```

Prepare the spatial part of a perturbed distribution function.

In general, this function returns

`` S = 1 + \\sum_m A_m \\cos ( k_m x ), ``

where the index ``m`` denotes the perturbed mode, `` k_m = 2\\pi m / L_x `` is the wavenumber of each mode,
and ``A_m`` is the corresponding amplitude.

# Notes
* This function accepts a single perturbed mode or a spectrum, but the number of perturbed modes
  must match with the number of amplitudes.
* The variable `modes` must be an `Integer`, a `Tuple` of `Integer`s or an `Array` of `Integer`s.
* The variable `amplitudes` must be an `Real`, a `Tuple` of `Real`s or an `Array` of `Real`s.
* In the case of using a 2-dimensional box, the keyword `dim` may be used to select which
  dimension should be perturbed.

# Examples
```jldoctest; setup = :(using Vlasova )
julia> box = Box(Nx = 256, Nv = 512, Lx = 5pi, vmin = -6, vmax = 6 );

julia> p = cosine_perturbation1d(box, modes = [1, 3], amplitudes = [1e-3, 1.5e-3] );

julia> sum(p) ≈ box.Nx[1]
true
```
"""
function cosine_perturbation1d(box::Box;
                               modes::Union{T, NTuple{N, T} where N, Array{T, 1}} where T <: Integer,
                               amplitudes::Union{T, NTuple{N, T} where N, Array{T, 1}} where T <: Real,
                               dim::Integer = 1)

    @assert length(modes) == length(amplitudes) "The perturbations and amplitudes should have the same length."
    @assert !in(0, modes) "The mode zero (0) can not be perturbed."

    x = box.x[dim]              # Space
    Lx = box.Lx[dim]            # Space length

    k = Float64[modes...] * 2pi / Lx # Perturbation spectrum
    A = Float64[amplitudes...]             # Perturbation amplitude(s)

    return 1 .+ reducedims(sum, A' .* cos.( k' .* x ), dims = 2)
end


"""
```julia
center(A::Array{T} where T <: Number; setmean = 0.0, dims = nothing)
```

Return an array consisting of `A` shifted so its new mean is `setmean`.
By default, the new mean value is `0.0` and the mean is taken along all dimensions
of `A`.

# Notes
* The returned array will always be of type `Float64` or `Complex{Float64}`.
* For an in-place version of this function, see [`center!`](@ref).
"""
function center(A::Array{T} where T <: Number; dims = nothing, setmean = 0.0)
    if isnothing(dims)
        M = Statistics.mean(A) - setmean
    else
        M = Statistics.mean(A, dims = dims) .- setmean
    end

    B = similar(A, promote_type( eltype(A), eltype(M) ) )

    @. B = A - M

    return B
end

"""
```julia
center!(A::Array{T} where T <: Union{Float64, Complex{Float64}}; dims = nothing, setmean = 0.0)
```
Shift the array `A` so that its new mean is `setmean`.
By default, the new mean value is `0.0` and the mean is taken along all dimensions of `A`.

# Notes
* The mean will always be of type `Float64` or `Complex{Float64}`. If the type of the elements
  of `A` is not one of either, the in-place operation will fail.
* For the not-in-place version of this function, see [`center`](@ref).

"""
function center!(A::Array{T} where T <: Union{Float64, Complex{Float64}}; dims = nothing, setmean = 0.0)
    if isnothing(dims)
        M = Statistics.mean(A) - setmean
    else
        M = Statistics.mean(A, dims = dims) .- setmean
    end

    @. A = A - M

    return nothing
end
