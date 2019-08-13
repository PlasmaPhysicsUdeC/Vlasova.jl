"""
    Perform a space advection efficiently by using the information stored in an
    object of type SpaceAdvection.
"""
function (space_advection::SpaceAdvection)(plasma::Plasma; advection_number::Integer = 1)

    for s in 1:plasma.number_of_species

        LinearAlgebra.mul!(space_advection.transformed_DF,
                           space_advection.plan, plasma.species[s].distribution)

        # Here happens the magic!
        _space_advection!( space_advection.transformed_DF,
                           space_advection.shift[advection_number][s] )

        LinearAlgebra.ldiv!(plasma.species[s].distribution,
                            space_advection.plan, space_advection.transformed_DF)
    end

    return 0;
end


@inline function _space_advection!(transformed_DF, shift)

    @. transformed_DF *= shift

    return 0;
end
