"""
Collects the total charge density of an element of type plasma
"""
function get_density(plasma::Plasma)
    
    chargedensity = zeros(Float64, plasma.box.Nx)
    for s in plasma.specie_axis
        chargedensity .+= plasma.species[s].charge * ( prod(plasma.box.dv) *
                                                       dropdims( sum( plasma.species[s].distribution,
                                                                      dims = plasma.box.velocity_dims ),
                                                                 dims = plasma.box.velocity_dims ) )
    end
    return chargedensity .- mean(chargedensity)
end

"""
Same as get_density but in place!
"""
function get_density!(chargedensity, plasma::Plasma)
    #=
    Same as get_density but in place
    =#
    
    chargedensity .= plasma.species[1].charge * ( prod(plasma.box.dv) *
                                                 dropdims( sum( plasma.species[1].distribution,
                                                                dims = plasma.box.velocity_dims ),
                                                           dims = plasma.box.velocity_dims ) )
    for s in 2:plasma.number_of_species
        chargedensity .+= plasma.species[s].charge * ( prod(plasma.box.dv) *
                                                       dropdims( sum( plasma.species[s].distribution,
                                                                      dims = plasma.box.velocity_dims ),
                                                                 dims = plasma.box.velocity_dims ) )
    end
    chargedensity .= chargedensity .- mean(chargedensity)
    return 0;
end


"""
Returns an array, where the i-th component corresponds to the kinetic energy of the
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
                                                  plasma.box,
                                                  temperature = plasma.species[i].temperature )
    end
    
    return kinetic_energies
end
