export notify, string2file, mean, suppress_stdout, hasnan
export @hasnan

# TODO: make a progressbar function
"""
    Prints a string both to the screen and to a file.
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
    Writes a string to a file

    Requires:
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
function mean(array::Array)::Number
    return sum(array) / length( array );
end

"""
    Executes a block of code without printing anything to screen

    Requires:
    * codeblock
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

#TODO check this
function halt_nan(var)
    hasnan = @hasnan var
    @assert !hasnan "Nan's found. Process halted."
end
