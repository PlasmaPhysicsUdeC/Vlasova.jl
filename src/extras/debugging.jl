"""
```julia
enable_debugging()
```

Enable debugging for the Vlasova module.

# Effects
* `TimerOutputs` (Package [`TimerOutputs.jl`](https://github.com/KristofferC/TimerOutputs.jl)):
  * Set a global `TimerOutput` called `timer` (unexported).
    Displaying it on the REPL (i.e. writing `Vlasova.timer`) will show the measured values.
  * Enable all the `@timeit_debug` annotations. This involves code recompilation, so it could take a while.

# Notes
* There is no option to disable debugging. This is accomplished by closing and re-starting the julia session.

"""
function enable_debugging()
    # Put everything debug-related in-here

    # Define timer and enable @timeit_debug
    global timer = TimerOutputs.TimerOutput()
    TimerOutputs.enable_debug_timings( Vlasova )

    return nothing
end
