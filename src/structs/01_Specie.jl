"""
    Defines a charged specie in a Vlasovian plasma.
    Contains and requires the following properties of a charged specie:
        * name :: String
        * charge :: Real
        * mass :: Real
        * temperature :: Real
        * distribution function :: Array{Real}
        where the distribution function is an array of 2N dimensions,
        and N is the number of spatial dimensions
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
    dim = Int( length(size(s.distribution)) / 2 )
    d = "
    ---
    $dim-dimensional Vlasova Specie.
    ---
    
    name = $(s.name)
    charge = $(s.charge)
    mass = $(s.mass)
    temperature = $(s.temperature)
    distribution size = $(size(s.distribution))"
    println(d)
end
