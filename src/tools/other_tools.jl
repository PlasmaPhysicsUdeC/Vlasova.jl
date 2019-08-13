export notify,
    string2file,
    mean,
    reducedims,
    suppress_stdout,
    hasnan,
    adiabatic_cutoff,
    outer,
    ⊗

export @hasnan

# TODO: make a progressbar function
"""
    Print a string both to the screen and to a file.
    ** If a global variable disable_notify = true is defined, notify will do nothing **

    Required:
    * string:: String to print

    Optional, keyword:
    * filename: String with the name of a file.
                If not specified, notify is equivalent to println()

    * mode: String specifying the mode under which filename will be open.
            [Options: "a" (append, default), "w" (write) ]
"""
function notify(string::String; filename::String = "/", mode::String = "a")
    (@isdefined disable_notify) && disable_notify ? (return nothing) : nothing

    println(string)
    (filename != "/") ? string2file(filename, string*"\n", mode) : nothing
    return nothing
end

"""
    Write a string to a file

    Required:
    * filename: String. The name of the target file
    * string: String to write to filename

    Optional
    * mode: Mode under which open filename.
            [Options: "r" (read), "w" (write), "a" (append, default), "r+" (special read & write)]
"""
function string2file(filename::String, string::String, mode::String = "a")
    open(filename, mode) do f
        write(f, string)
    end
    return nothing
end

"""
    Get the mean value of an N-dimensional array

    Requires:
    * array: Array

    Returns
    * mean: Float64
"""
function mean(array::Array; dims = (0))
    if dims == (0)
        return sum(array) / prod(size( array ));
    else
        return sum(array, dims = dims) / prod(size( array )[ [dims...] ]);
    end
end

# Change reducedims for a @squeezing macro
"""
    Performs a reduction function, f,  over an array, A, along the dimensions, dims,
    and drop the dimensions reduced.

    In general,
    ```
    reducedims( f, A, dims = dims)
    ```
    is equivalent to
    ```
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
    Execute a block of code without printing anything to screen
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
    Test whether some element of var (or var itself) is a NaN
"""
macro hasnan(var)
    quote
        findfirst(isnan.($var)) != nothing
    end
end


"""
    Test whether some element of var (or var itself) is a NaN
"""
function hasnan(var)
    return findfirst(isnan.(var)) != nothing
end

"""
    Return 1.0 for times before cutoff_time, adiabatically go to 0.0 until cutoff_time + cutoff_delay and return 0.0 afterwards

    The shape used to go from 1.0 to 0.0 is
    ```math
    f(x) =  \\frac{( 1 + cos(x) )}{2}
    ```
    with ```x \\in [0, pi]```.
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

raw"""
    Generalized external operation between two arrays, A and B.
    
    In general, `outer(A, B, f)` returns `R` such that
    ```math
    R_{ij} = f(A_i, B_j)
    ```

    By default, `f` is the usual multiplication, so that `outer` is the external product.

    ```@example
    A = [1.0, 2.0, 3.0, 4.0];
    B = [0.0, 1.0];
    outer(A, B, *)
    outer(A, B, -)
    ```

    # Note
    This function is not optimized to be fast, but only to offer a comfortable syntax.
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
Infix form of the external product between two arrays

See ["outer"](@ref)
"""
⊗ = outer
