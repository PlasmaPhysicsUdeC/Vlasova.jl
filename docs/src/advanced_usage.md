## Injecting code to the integrator

Sometimes one may need Vlasova to compute and save some quantities apart from the common ones. Since it is very difficult to provide options for all the possible situations, a function was implemented to allow the final user to directly inject code into the main integrator. However, this function is a double sided sword, and it should be handled as such.

```@docs
Vlasova.inject_to_integrator
```

<!-- # TODO: Example of injecting code to save something -->

## Debugging
To perform basic debugging measures, Vlasova may use [`TimerOutputs.jl`](https://github.com/KristofferC/TimerOutputs.jl) to analyze the time and memory allocations of each part of the main integrator.

By default, all debugging measurements are disabled yielding zero overhead, but they may be activated by calling the unexported function

```@docs
Vlasova.enable_debugging()
```

which sets a `TimerOutput` called `Vlasova.timer`, collecting the debugging information as the simulation runs. To view this information, the timer must be displayed in the REPL by simply calling

```julia
julia> Vlasova.timer
```
