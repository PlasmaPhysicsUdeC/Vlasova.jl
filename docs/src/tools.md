## Basic Structs
```@autodocs
Modules = [Vlasova]
Pages = ["01_Box.jl", "01_Specie.jl", "02_Plasma.jl"]
Private = false
Public = true
```

## Symplectic Integrators

# TODO: Cite verlet integrators
Vlasova implements a way to define arbitrary gradient and non-gradient symplectic integrators with the function [`VlasovaIntegrator`](@ref). Nonetheless, a few integrators come defined by default, and those integrators are the following (unless explicitly cited, they were taken from [`Omelyan (2003)`](https://www.sciencedirect.com/science/article/pii/S0010465502007543) ).

```@autodocs
Modules = [Vlasova]
Pages = ["integrators.jl"]
Private = false
Public = true
```

## Distributions
```@autodocs
Modules = [Vlasova]
Pages = ["distributions.jl"]
Private = false
Public = true
```

## Field Analysis
<!-- Add simulation/post-simulation tools-->

## Mathematical operations

## Other Tools
# TODO: Separate the mathematical operator from here, and put them in their section

```@autodocs
Modules = [Vlasova]
Pages = ["other_tools.jl"]
Private = false
Public = true
```
