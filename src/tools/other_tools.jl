export notify, string2file, mean, suppress_stdout

# TODO: make a progressbar function

function notify(string::String; fileName::String = "/", mode::String = "a")
    #=
    Prints string both to the screen and to a file.

    * If no name is provided for the file, notify() is equivalent to println()
    * If disable_notify = true is defined, notify() does nothing.
    =#
    (@isdefined disable_notify) && disable_notify ? (return nothing) : nothing

    println(string)
    (fileName != "/") ? string2file(fileName, string*"\n", mode) : nothing
    return nothing
end

function string2file(fileName::String, string::String, mode::String = "a")
    #=
    Appends "string" to the end of a file "fileName" and immediately closes it.
    Requires:
    * fileName
    * string

    [optional]
    * mode: "r", read; "w", write; "a", append; "r+", special read & write
    =#
    open(fileName, mode) do f
        write(f, string)
    end
    return nothing
end


function mean(array::Array)::Number
    #=
    Get the mean valuie of an N-dimensional array

    Requires:
    * array

    =#
    return sum(array) / length( array );
end


macro suppress_stdout(codeBlock)
    #=
    Executes a block of code without printing anything to screen

    Requires:
    * codeBlock
    =#
    quote
        if ccall(:jl_generating_output, Cint, ()) == 0
            original_stdout = stdout
            out_rd, out_wr = redirect_stdout()
            out_reader = @async read(out_rd, String)
        end

        try
            $(esc(codeBlock))
        finally
            if ccall(:jl_generating_output, Cint, ()) == 0
                redirect_stdout(original_stdout)
                close(out_wr)
            end
        end
    end
end
