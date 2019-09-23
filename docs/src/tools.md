## Basic Structs
```@autodocs
Modules = [Vlasova]
Pages = ["01_Box.jl", "01_Specie.jl", "02_Plasma.jl"]
Private = false
Public = true
```

## Symplectic Integrators

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

## Result analysis
Since the electrostatic regime is assumed, it is generally easier to work with the charge density, instead of the electric field. For this reason Vlasova saves the charge density of the simulations, and consequently, all the tools to analyze the fields take the charge density as input.

```@autodocs
Modules = [Vlasova]
Pages = ["result_analysis.jl"]
Private = false
Public = true
```

## Convenience and mathematical tools
```@autodocs
Modules = [Vlasova]
Pages = ["convenience_tools.jl"]
Private = false
Public = true
```

```@autodocs
Modules = [Vlasova]
Pages = ["fourier_tools.jl"]
Private = false
Public = true
```

## Other Tools

```@autodocs
Modules = [Vlasova]
Pages = ["other_tools.jl"]
Private = false
Public = true
```
