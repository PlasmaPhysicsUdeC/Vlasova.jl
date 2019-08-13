"""
Collect the total charge density of an element of type plasma.
Prefer in-place version.
"""
function get_density(plasma::Plasma)

    # Allocate array
    chargedensity = Array{Float64}(undef, plasma.box.Nx)

    # Call in-place version
    get_density!(chargedensity, plasma)

    return chargedensity
end

"""
Same as get_density but in place!
"""
function get_density!(chargedensity, plasma::Plasma)
    # Replace chargedensity content
    chargedensity .= plasma.species[1].charge * ( prod(plasma.box.dv) *
                                                 dropdims( sum( plasma.species[1].distribution,
                                                                dims = plasma.box.velocity_dims ),
                                                           dims = plasma.box.velocity_dims ) )
    # Loop over species
    for s in 2:plasma.number_of_species
        chargedensity .+= plasma.species[s].charge * ( prod(plasma.box.dv) *
                                                       dropdims( sum( plasma.species[s].distribution,
                                                                      dims = plasma.box.velocity_dims ),
                                                                 dims = plasma.box.velocity_dims ) )
    end

    # Ensure quasi-neutrality
    chargedensity .= chargedensity .- mean(chargedensity)
    return 0;
end


"""
Return an array of kinetic energies, where the i-th component corresponds to the i-th charged specie of the plasma.
"""
function get_kinetic_energies(plasma::Plasma)

    return  [ get_kinetic_energy( plasma.species[i].distribution,
                                  plasma.box,
                                  temperature = plasma.species[i].temperature ) for i in plasma.specie_axis ]
end
