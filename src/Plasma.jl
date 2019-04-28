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
mutable struct Specie
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
           ) = new( name,
                    Float64(charge),
                    Float64(mass),
                    Float64(temperature),
                    Float64.(distribution) )
end

"""
    Defines a Vlasovian Plasma, and contains:
        * species: An array of Specie elements
        * box: The parameters of the simulation, of type Box
        * number_of_species: The number of charged species contained in the plasma
        * specie_axis: A tuple of indices, to iterate over each specie.
        
        To be defined, it required to specify species and box:
            julia> parameters = Box( ... );
            julia> plasma = Plasma([Specie(...),
                                    Specie(...)],
                                    parameters );
        *** This struct is dependent on the Specie and Box structs ***
        """
mutable struct Plasma
    species::Array{Specie}
    box::Box
    number_of_species::Int64
    specie_axis::Base.OneTo{Int64}
    
    # Constructors
    # Calculate derived quantities from the array of species
    Plasma(s::Array{Specie},
           p::Box ) = new(s,
                          p,
                          size(s, 1),
                          Base.OneTo(size(s, 1))
                          )
end

"""
Collects the total charge density of an element of type plasma
"""
function get_density(plasma::Plasma)
    
    chargedensity = zeros(Float64, plasma.box.Nx)
    for s in plasma.specie_axis
        chargedensity .+= plasma.species[s].charge * ( prod(plasma.box.dv) *
                                                        dropdims( sum( plasma.specie[s].distribution,
                                                                       dims = plasma.box.space_dims ),
                                                                  dims = plasma.box.space_dims ) )
    end
    return chargedensity
end

"""
Same as get_density but in place!
"""
function get_density!(chargedensity, plasma::Plasma)
    #=
    Same as get_density but in place
    =#
    
    chargedensity = plasma.species[1].charge * ( prod(plasma.box.dv) *
                                                 dropdims( sum( plasma.species[1].distribution,
                                                                dims = plasma.box.space_dims ),
                                                           dims = plasma.box.space_dims ) )
    
    for s in 2:plasma.number_of_species
        chargedensity .+= plasma.species[s].charge * ( prod(plasma.box.dv) *
                                                       dropdims( sum( plasma.species[s].distribution,
                                                                      dims = plasma.box.space_dims ),
                                                                 dims = plasma.box.space_dims ) )
    end
    return 0;
end


"""
Returns an array, where the i-th component corresponsd to the kinetic energy of the
i-th charged specie of the plasma.
"""
function get_kinetic_energies(plasma::Plasma)
    #=
    Returns the respective kinetic energies of each specie in an array
    Requires:
    * Species
    =#

    kinetic_energies = Array{Float64}(undef, plasma.number_of_species)

    for i in plasma.specie_axis
        kinetic_energies[i] = get_kinetic_energy( plasma.species[i].distribution,
                                                  plasma.box.v,
                                                  plasma.box.dx,
                                                  temperature = plasma.species[i].temperature )
    end
    
    return kinetic_energies
end
