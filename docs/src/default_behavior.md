## Default behavior

This page lists all the variables that may be customized to control the behavior of [`integrate!`](@ref).

This variables are intentionally unexported, meaning that they must be accessed as `Vlasova.variable_name`, and they can only be changed from the inside of the module defined by Vlasova. To accomplish the latter, the macro [`@vlasova`](@ref) has been implemented.

```@autodocs
Modules = [Vlasova]
Pages = ["defaults.jl"]
Private = true
```
