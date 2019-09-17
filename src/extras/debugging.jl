"""
```julia
enable_debugging()
```

Enable debugging for the Vlasova module.

# Notes
* This function should recompile all code involving debugging annotations. For this reason it may take a while.

"""
function enable_debugging()
    # Put everything debug-related in-here

    # Enable @timeit_debug
    TimerOutputs.enable_debug_timings( Vlasova )
end
