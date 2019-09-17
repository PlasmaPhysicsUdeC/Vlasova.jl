"""
```julia
enable_debugging()
```

Enable debugging for the Vlasova module.

# Effects
* `TimerOutputs`:
  * Set a global `TimerOutput` called `timer`. Calling in on the REPL shows the measured values.
  * Enable all the `@timeit_debug` annotations. This involves code recompilation, so it could take a while.

"""
function enable_debugging()
    # Put everything debug-related in-here

    # Define timer and enable @timeit_debug
    global timer = TimerOutputs.TimerOutput()
    TimerOutputs.enable_debug_timings( Vlasova )
end
