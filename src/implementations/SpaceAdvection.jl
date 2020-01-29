"""
    Perform a space advection efficiently by using the information stored in an
    object of type SpaceAdvection.
"""
function (space_advection::SpaceAdvection)(plasma::Plasma; advection_number::Integer = 1)

    for s in 1:plasma.number_of_species

        TimerOutputs.@timeit_debug timer "Fourier Transform" begin
            LinearAlgebra.mul!(space_advection.transformed_DF,
                               space_advection.plan, plasma.species[s].distribution)
        end

        # Here happens the magic!
        TimerOutputs.@timeit_debug timer "Apply propagator" begin
            _space_advection!( space_advection.transformed_DF,
                               space_advection.shift[advection_number][s] )
        end

        TimerOutputs.@timeit_debug timer "inverse Fourier transform" begin
            LinearAlgebra.ldiv!(plasma.species[s].distribution,
                                space_advection.plan, space_advection.transformed_DF)
        end
    end

    return nothing;
end


@inline function _space_advection!(transformed_DF, shift)

    Strided.@strided transformed_DF .*= shift

    return nothing;
end
