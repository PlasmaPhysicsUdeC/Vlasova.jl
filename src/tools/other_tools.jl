export @vlasova,
    @hasnan,
    @suppress_stdout,
    hasnan,
    adiabatic_cutoff,
    sim_from_file

"""
```julia
@vlasova codeblock
```

Evaluate `codeblock` inside the scope of the Vlasova module.

# Examples

To change the default value of an inner variable of Vlasova, such as `continue_from_backup`, use

```jldoctest; setup = :(using Vlasova)
julia> @vlasova continue_from_backup = true
2
```
"""
macro vlasova( codeblock )
    Vlasova.eval( codeblock )
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
