"""
```julia
Specie(;name::String, charge::Real, mass:Real, temperature::Real, distribution::Array{T} where T :< Real)
```

Create a container representing one charged species in a Vlasova plasma.

# Notes
* All `Real`s will be converted to `Float64` internally.
"""
struct Specie                   # TODO: Accept thermal_velocity or temperature. (rewrite integrator to use vth)
    name::String
    charge::Float64
    mass::Float64
    temperature::Float64
    distribution::Array{Float64}

    # Constructors
    # Take reals, but turn them into Float64
    _Specie(name::String,
            charge::Real,
            mass::Real,
            temperature::Real,
            distribution::Array{T} where T <: Real
            ) = begin
                dist_dim = length( size(distribution) )
                @assert dist_dim > 1 "The distribution provided is 1-dimensional. It should have 2n dimensions where n >= 1"
                @assert mod(dist_dim, 2) == 0 "The distribution provided should be 2n dimensional (n dimensions in space and n dimensions in velocity)"
                new( name,
                     Float64(charge),
                     Float64(mass),
                     Float64(temperature),
                     Float64.(distribution) )
            end

    Specie(args...
           ;name::String,
           charge::Real,
           mass::Real,
           temperature::Real,
           distribution::Array{T} where T <: Real) = begin
               @assert args == () "Specie should not be called using non-keyword arguments."
               _Specie(name,
                       charge,
                       mass,
                       temperature,
                       distribution)
           end
end


function Base.display(s::Specie)
    dim = div(length(size(s.distribution)), 2)
    d = "
    ---
    $dim-dimensional Vlasova Specie.
    ---

    name = $(s.name)
    charge = $(s.charge)
    mass = $(s.mass)
    temperature = $(s.temperature)
    distribution size = $(size(s.distribution))
    "
    print(d)
end

function Base.display(s::Array{Specie})
    dim = div(length(size(s[1].distribution)), 2)
    sp_list = ""
    for i in 1:size(s, 1)
        sp_list *= "$i) "*s[i].name*"\n    "
    end
    d = "
    ---
    Array of $dim-dimensional Vlasova Species.
    ---

    Species contained:
    $sp_list "
    print(d)
end
