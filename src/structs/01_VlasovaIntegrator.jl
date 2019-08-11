struct VlasovaIntegrator
    sequence::String
    coefficients::Array{T} where T <: AbstractFloat
    gradient_coefficients::Array{T} where T <: AbstractFloat
    merge_last_advection::Bool # TODO: MERGE

    VlasovaIntegrator( sequence::String,
                       coefficients::Array{T} where T <: AbstractFloat;
                       gradient_coefficients::Array{T} where T <: AbstractFloat = Float64[]
                       ) = begin
                           @assert length(sequence) == length(coefficients) "The number of coefficients do not correspond to the number of advections"

                           for i in sequence
                               @assert (i in "ABC") "The integrator sequence must only contain A's, B's and C's"
                           end
                           
                           C_locations = [i == 'C' for i in sequence ]
                           @assert sum(C_locations) == length(gradient_coefficients) "The number of gradient steps (C's) and the number of gradient coefficients do not match"
                           
                           pos_locations = [i == 'A' for i in sequence ]
                           pos_ind = findall( pos_locations )
                           vel_ind = findall( .!pos_locations )
                           
                           @assert isapprox(sum(coefficients[pos_ind]), 1) "The position coefficients do not add up to 1"
                           @assert isapprox(sum(coefficients[vel_ind]), 1) "The velocity coefficients do not add up to 1"
                           
                           for i in 2:length( sequence )
                               same = sequence[i-1] == sequence[i]
                               @assert !same "Advections of kind A and B/C must be interleaved."
                           end
                           @assert !any(occursin.(["CB", "BC"], sequence)) "There can't be succesive B's and C's in the integrator sequence" 
                           # if merge_last_advection
                           #     @assert (sequence[1] == sequence[end]) "Can't merge first and last advections if they are not of the same kind"
                           # end
                           new(sequence,
                               coefficients,
                               gradient_coefficients,
                               false)
                       end
    
    VlasovaIntegrator(sequence::String,
                      coefficients::NTuple{N, T} where N where T <: AbstractFloat;
                      gradient_coefficients = Float64[] ) =
                          VlasovaIntegrator(
                              sequence,
                              Array([coefficients...]),
                              gradient_coefficients = gradient_coefficients
                          )
end
