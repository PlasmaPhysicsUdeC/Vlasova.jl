export @vlasova,
    @hasnan,
    @suppress_stdout,
    hasnan,
    reducedims,
    adiabatic_cutoff,
    outer,
    ⊗,
    cosine_perturbation1d,
    sim_from_file,
    find_local_maxima,
    find_exponential_growth

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
julia> @suppress_stdout println(2) # Prints nothing

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

Test whether some element of `var` (or itself) is a `NaN`.

# Notes
* There is also a function version [`hasnan`](@ref).

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
    return quote
        findfirst(isnan.($var)) != nothing
    end
end

"""
```julia
hasnan( var )
```

Test whether some element of `var` (or itself) is a `NaN`.

# Notes
* There is also a macro version [`@hasnan`](@ref).

# Examples
```jldoctest; setup = :(using Vlasova)
julia> hasnan( NaN )
true

julia> hasnan( [1.0, NaN] )
true

julia> hasnan( [1, 2, Inf] )
false
```
"""
function hasnan(var)
    return findfirst(isnan.(var)) != nothing
end

"""
```julia
adiabatic_cutoff(time::Real; cutoff_time::Real, cutoff_delay::Real)
```

Return ``1.0`` if  `time < cutoff_time`, adiabatically go to ``0.0``
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
sim_from_file(parfile::String)
```

Start a simulation from a file that defines at least the `plasma`, `final_time`, and `dt` (with those names).

# Notes

* When the variable [`Vlasova.data_path`](@ref) is set, this function will make a copy of the file
  of parameters into that path.
* This function was developed to simplify the workflow by just having a file defining tha basic parameters
  and running the simulation in one instruction, while having the possibility to save the results and the
  file with parameters itself.


# Examples

Supposing that there is a file "parameters.jl" in the current path, and it defines the variables
specified, a simulation may be performed as it follows

```julia
julia> using Vlasova

julia> sim_from_file("parameters.jl")
Preparing integrator. This may take a while...
Starting integration @ (...)
Progress: 100%[======================================] Time: (...)
```

The same as before could be performed in the following one-liner directly from the terminal

```bash
\$ julia -e 'using Vlasova; sim_from_file("parameters.jl\")'
```
which allows the use of tools such as nohup.
"""
function sim_from_file(parfile::String)
    # Check existence of the required variables
    for var in [:plasma, :final_time, :dt]
        @assert (@isdefined var) "The parameters file must define the $(String(var))."
    end

    # Evaluate in the Main module
    @eval Main begin
        # Include parameters
        include( $parfile )
        # Copy parameters to folder file if defined
        mkpath( Vlasova.data_path )
        !isempty( Vlasova.data_path ) ? cp( $parfile, joinpath( Vlasova.data_path, $parfile), force = true) : nothing
        # Start integrator
        integrate!(plasma, final_time, dt)
    end
end


"""
```julia
find_local_maxima(a::AbstractArray{T, 1} where T <: Real; discard_borders = false )
```

Find the indices of all the local maxima of a real-valued one-dimensional `Array`.

# Notes
* If `discard_borders = true`, the values `a[1]` and `a[end]` are not tested to be local maxima.
  This option is useful when the local maxima are wanted to fit the data.

# Examples

```jldoctest; setup = :(using Vlasova)
julia> x = range(0.0, stop = 15pi, step = 0.1);

julia> y = cos.(x); # Has a maximum on the first index.

julia> maxima = find_local_maxima( y );

julia> 1 in maxima  # Is 1 on the maxima?
true

julia> maxima = find_local_maxima( y, discard_borders = true);

julia> 1 in maxima # Is 1 on the maxima?
false

"""
function find_local_maxima(a::AbstractArray{T, 1} where T <: Real; discard_borders = false  )
    sa = size(a, 1)
    maxinds = Int64[]
    if sa <= 1
        discard_borders ? nothing : ( maxinds = copy(a) )
    else
        # Middle indices
        for i in 2:(sa-1)
            a[i-1] < a[i] > a[i+1] ? push!(maxinds, i) : nothing
        end

        if !discard_borders
            # First index
            a[1] > a[2] ? pushfirst!(maxinds, 1) : nothing
            # End index
            a[end] > a[end-1] ? push!(maxinds, sa) : nothing
        end
    end

    return maxinds
end


"""
```julia
find_exponential_growth( x, y; interval = [-Inf, Inf] )
```

Fit the function

``y(x) = A \\exp ( B x  )``

 over the interval given to the `Array`s `x` and `y`, and return the growth rate, ``B``.

# Notes
* `x` and `y` must be `Array`s of the same length.
* The exponential growth returned is the additive inverse of the damping rate.

# Examples
```jldoctest; setup = :(using Vlasova)
julia> x = range(0.0, stop = 5pi, step = 0.1);

julia> A, B = rand(2);

julia> y = A * exp.(B * x); # Fake data

julia> growth = find_exponential_growth(x, y);

julia> growth ≈ B
true

```
"""
function find_exponential_growth( x, y; interval = [-Inf, Inf] )
    idx = findall( interval[1] .<= x .<= interval[2] )
    return CurveFit.exp_fit(x[idx], y[idx])[2]
end
