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
struct Specie
    name::String
    charge::Float64
    mass::Float64
    temperature::Float64
    distribution::Array{Float64}

    # Constructors
    # Take reals, but turn them into Float64
    Specie(name::String,
           charge::Real,
           mass::Real,
           temperature::Real,
           distribution::Array{T} where T <: Real
           ) = begin
               @assert length( size(distribution) ) > 1 "The distribution provided is 1-dimensional. It should have 2n dimensions where n >= 1"
               new( name,
                    Float64(charge),
                    Float64(mass),
                    Float64(temperature),
                    Float64.(distribution) )
           end
end
