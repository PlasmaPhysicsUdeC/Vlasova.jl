"""
Executed at every iteration of the integrator.
"""
function (d::DataSaver)(plasma::Plasma, t::Integer )
    if d.save_data
        # Saved every iteration
        it = t - d.last_iteration_saved
        d.chargedensity[plasma.box.space_axes..., it] = get_density( plasma )
        d.kinetic_energy[it, :] = [ get_kinetic_energy(plasma.box, plasma.species[s])
                                    for s in plasma.specie_axis ]

        # Save distribution ?
        ( t in d.save_distribution_axis ) ? save_distribution(d, plasma) : nothing

        # is checkpoint ?
        if t in d.checkpoint_axis
            save_to_disk(d, plasma, t)
        end
    end

    return nothing;
end

"""
Flush the values accumulated on a `DataSaver` object to disk.
"""
function save_to_disk(d::DataSaver, p::Plasma, t::Integer)
    last_it = d.last_iteration_saved

    # For the common files
    fid = HDF5.h5open(joinpath(d.path, "shared_data.h5"), "r+")
    fid["chargedensity"][p.box.space_axes..., (last_it+1):t] = d.chargedensity[p.box.space_axes..., 1:(t-last_it)]
    fid["total_kinetic_energy"][(last_it + 1):t] = reducedims(sum,
                                                              d.kinetic_energy[1:(t-last_it), :],
                                                              dims = 2 )
    HDF5.close(fid)

    for s in 1:p.number_of_species
        fid = HDF5.h5open(joinpath(d.path, p.species[s].name*".h5"), "r+")
        s == 1 ? (fid["last_iteration_saved"][1] = t) : nothing
        fid["distribution"][p.box.distribution_axes..., length(d.save_distribution_axis)] = p.species[s].distribution
        HDF5.close(fid)
    end

    # Update indices
    d.last_iteration_saved = t
    d.checkpoints_reached += 1

    # At the last iteration, specify that the files are complete
    if t == d.checkpoint_axis[end]
        fid = HDF5.h5open(joinpath(d.path, "shared_data.h5"), "r+")
        # TODO: There is a Bug by now on HDF5 when saving strings: https://github.com/JuliaLang/julia/issues/32752
        HDF5.attrs(fid)["Completed"] = true#"Complete"
        HDF5.close(fid)
    end

    return nothing;
end

function save_distribution(d::DataSaver, p::Plasma)
    for s in 1:p.number_of_species
        fid = HDF5.h5open( joinpath(d.path, p.species[s].name*".h5"), "r+")
        fid["distribution"][p.box.distribution_axes..., d.last_distribution_saved + 1 ] = p.species[s].distribution
        HDF5.close(fid)
    end
    d.last_distribution_saved += 1

    return nothing;
end
