export @vlasova,
    @hasnan,
    @suppress_stdout,
    reducedims,
    hasnan,
    adiabatic_cutoff,
    outer,
    ⊗

"""
```julia
@vlasova codeblock
```

Evaluate `codeblock` inside the scope of the Vlasova module.

# Examples

To change the default value of an inner variable of Vlasova, such as `NUM_THREADS`, use

```jldoctest; setup = :(using Vlasova)
julia> @vlasova NUM_THREADS = 2
2
```
"""
macro vlasova( codeblock )
    Vlasova.eval( codeblock )
end

"""
```julia
reducedims(f::Function, A::Array; dims = (0))
```

Perform a reduction function, `f`,  over an array, `A`, along the dimension(s), `dims`,
and drop the dimensions reduced.

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
@supress_stdout codeblock
```

Evaluate `codeblock` without printing anything to screen.

# Examples
```jldoctest; setup = :(using Vlasova)
julia> @suppress_stdout println(2)

```
"""
macro suppress_stdout(codeblock)
    quote
        if ccall(:jl_generating_output, Cint, ()) == 0
            original_stdout = stdout
            out_rd, out_wr = redirect_stdout()
            out_reader = @async read(out_rd, String)
        end

        try
            $(esc(codeblock))
        finally
            if ccall(:jl_generating_output, Cint, ()) == 0
                redirect_stdout(original_stdout)
                close(out_wr)
            end
        end
    end
end

"""
```julia
@hasnan var
```

Test whether some element of `var` (or itself) is a `NaN`

# Examples
```jldoctest; setup = :(using Vlasova)
julia> @hasnan NaN
true

julia> @hasnan [1.0, NaN]
true

julia> @hasnan [1, 2, Inf]
false
```
"""
macro hasnan(var)
    quote
        findfirst(isnan.($var)) != nothing
    end
end


"""
Test whether some element of `var` (or itself) is a `NaN`
"""
function hasnan(var)
    return findfirst(isnan.(var)) != nothing
end

"""
```julia
adiabatic_cutoff(time::Real; cutoff_time::Real, cutoff_delay::Real)
```

Return ``1.0`` if  `time`<`cutoff_time`, adiabatically go to ``0.0``
 until `time = cutoff_time + cutoff_delay` and return `0.0` afterwards.

The shape used to go from ``1.0`` to ``0.0`` is

``f(x) =  \\frac{( 1 + \\cos(x) )}{2},`` with ``x \\in [0, \\pi]``.
"""
function adiabatic_cutoff(time::Real; cutoff_time::Real, cutoff_delay::Real)
    damp_time = time - cutoff_time
    if damp_time < 0
        return 1.0
    elseif damp_time < cutoff_delay
        return ( 1.0 + cos( pi * damp_time / cutoff_delay ) ) / 2
    else
        return 0.0
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
