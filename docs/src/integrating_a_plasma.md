## Integrating a plasma

At this point, it is assumed that a [`Plasma`](@ref) was created from a [`Box`](@ref) and an array of [`Specie`](@ref)s. As it was said earlier, the plasma may be evolved in time just by typing

```julia
integrate!(plasma, final_time, dt),
```

which will evolve `plasma` until time `final_time` using steps `dt`.

However, a lot more may be customized.

Currently, the variables that may be customized to control the behavior of [`integrate!`](@ref) are

## Default behavior
```@autodocs
Modules = [Vlasova]
Pages = ["defaults.jl"]
Private = true
```
