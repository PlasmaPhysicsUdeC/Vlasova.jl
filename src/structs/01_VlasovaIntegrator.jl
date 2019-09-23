"""
```julia
VlasovaIntegrator( sequence::String,
                   coefficients::Array{T} where T <: AbstractFloat;
                   gradient_coefficients::Array{T} where T <: AbstractFloat = Float64[]
                  )
```

Define a symplectic integrator to be used by Vlasova.

The convention used by this struct is the following:

* `sequence`: It is a string which contains as many characters as the number of advections performed
   for each time step. The characters may only be `'A'`, `'B'` or `'C'`, where:
   * `'A'`: Corresponds to a space advection.
   * `'B'`: Corresponds to a velocity advection.
   * `'C'`: Corresponds to a velocity advection with gradient correction.

* `coefficients`: Is an array of `Float64`s, where the \$ n\$-th element of the array is the
   advection coefficient of the \$n\$-th advection in `sequence`.

If the integrator has gradient force corrections, the keyword argument `gradient_coefficients` must be provided,
being an array where the \$n\$-th element is the coefficient of the force correction of the \$n\$-th gradient-corrected advection.

# Examples

The velocity form of the Verlet integrator [1], where the integration is split in
* a velocity advection for half a time step,
* a space advection for a whole time step, and
* a velocity advection for another half of a time step,

is defined as
```julia
verlet_velocity = VlasovaIntegrator("BAB", [0.5, 1.0, 0.5] )
```

In a similar manner, the fourth order integrator Chin A [2], is defined as
```julia
chin_A = Vlasovaintegrator("BACAB", [ 1/6, 3/8, 1/3, 1/4, 1/3, 3/8, 1/6 ],
                            gradient_coefficients = [ 1/48 ] )
```

# References
[1]: https://aip.scitation.org/doi/10.1063/1.442716
[2]: https://www.sciencedirect.com/science/article/pii/S0375960197000030
"""
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
