struct VlasovaIntegrator
    sequence::String
    coefficients::Array{T} where T <: AbstractFloat
    merge_last_advection::Bool # TODO: MERGE

    VlasovaIntegrator( sequence::String,
                       coefficients::Array{T} where T <: AbstractFloat
                       ) = begin
                           @assert length(sequence) == length(coefficients) "The number of coefficients do not correspond to the number of advections"

                           for i in sequence
                               @assert (i in "AB") "The integrator sequence must only contain A's and B's"
                           end

                           pos_locations = [i == 'A' for i in sequence ]
                           pos_ind = findall( pos_locations )
                           vel_ind = findall( .!pos_locations )
                           
                           @assert isapprox(sum(coefficients[pos_ind]), 1) "The position coefficients do not add up to 1"
                           @assert isapprox(sum(coefficients[vel_ind]), 1) "The velocity coefficients do not add up to 1"
                           
                           isA = sequence[1] == 'A'
                           for i in sequence[2:end]
                               old = isA
                               isA = i == 'A'
                               @assert xor(old, isA) "Advections of kind A and B must be interleaved."
                           end
                           # if merge_last_advection
                           #     @assert (sequence[1] == sequence[end]) "Can't merge first and last advections if they are not of the same kind"
                           # end
                           new(sequence,
                               coefficients,
                               false)
                       end
    VlasovaIntegrator(sequence::String,
                      coefficients::NTuple{N, T} where N where T <: AbstractFloat ) =
                          new(
                              sequence,
                              Array([coefficients...])
                              )
end
