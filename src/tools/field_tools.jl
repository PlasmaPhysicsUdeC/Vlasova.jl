export get_kinetic_energy,
    get_density,
    get_density!,
    get_electric_field,
    get_electrostatic_energy,
    get_power_per_mode,
    get_dispersion_relation

"""
```julia
get_density(box::Box, species::Specie)
```

Collect the charge density of a [`Specie`](@ref).
"""
@inline function get_density(box::Box, species::Specie)
    return species.charge * prod(box.dv) * reducedims(sum,
                                                      species.distribution,
                                                      dims = box.velocity_dims)
end

# TODO: Verify @ref pointing with specialization.
"""
```julia
get_density(plasma::Plasma)
```

Collect and return the total charge density of a [`Plasma`](@ref).
"""
function get_density(plasma::Plasma)

    # Allocate array
    chargedensity = Array{Float64}(undef, plasma.box.Nx)

    # Call in-place version
    get_density!(chargedensity, plasma)

    # Ensure quasi-neutrality
    chargedensity .= chargedensity .- Statistics.mean(chargedensity)

    return chargedensity
end

"""
```julia
get_density!(chargedensity::Array{Float64}, plasma::Plasma)
```

Collect the total `chargedensity` of a [`Plasma`](@ref) in place.
"""
function get_density!(chargedensity::Array{Float64}, plasma::Plasma)
    # Replace chargedensity content
    chargedensity .= get_density(plasma.box, plasma.species[1])

    # Loop over species
    for s in 2:plasma.number_of_species
        chargedensity .+= get_density(plasma.box, plasma.species[s])
    end

    # Ensure quasi-neutrality
    chargedensity .= chargedensity .- Statistics.mean(chargedensity)
    return nothing;
end

"""
```julia
get_kinetic_energy(plasma::Plasma)
```

Obtain the total kinetic energy of a [`Plasma`](@ref).
"""
function get_kinetic_energy(plasma::Plasma)
    kinen = 0.0
    for s in plasma.box.specie_axis
        kinen += get_kinetic_energy(plasma.box, plasma.species[s])
    end
    return kinen
end

"""
```julia
get_kinetic_energy(box::Box, species::Specie);
```

Obtain the kinetic energy from a [`Specie`](@ref) and a [`Box`](@ref).

"""
function get_kinetic_energy(box::Box, species::Specie);

    fred = reducedims(sum,
                      species.distribution,
                      dims = box.space_dims )

    s = 0.0
    for i in CartesianIndices( box.Nv ) # This reduction could be parallelized
        v2 = 0.0
        for d in box.dim_axis
            v2 += box.v[d][i[d]]^2
        end
        s += fred[i] * v2
    end

    return species.temperature * prod(box.dx) * prod(box.dv) * s / 2;
end

"""
```julia
get_electric_field(box::Box, chargedensity::Array{Float64})
```

Obtain the electric field from a charge distribution.

# Notes
* This function follows the usual Vlasova convention that scalars such
  as the `chargedensity`are `Array{Float64}` and the vectors such as
  the `electric_field` are `Array{Array{Float64}}`, where each component
  corresponds to the same component of the vector which is a scalar (`Array{Float64}`).
  This is true even in the 1-dimensional case.

* The `chargedensity` may depend upon time over the last axis, in which case
  each component of the `electric_field` will have the same dependence on time.

# Examples

```jldoctest; setup = :(using Vlasova)
julia> box = Box(Nx = 256, Nv = 512, Lx = 2pi, vmin = -6,  vmax = 6 );

julia> chargedensity = sin.( box.x[1] ); # Fake data FTW!

julia> E = get_electric_field(box, chargedensity); # Should be ``-cos(x)``

julia> E[1] ≈ -cos.( box.x[1] ) # E[1] since E is an Array{Array{Float64}}
true
```
"""
function get_electric_field(box::Box, chargedensity::Array{Float64})

    Nx2p1, fourier_axis = get_rfft_dims( box.x )

    k = rfft_wavevector( box.x )
    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0

    integrate = Array{Array{Complex{Float64}}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        integrate[d] = -1im ./ k2
        for i in fourier_axis
            integrate[d][ i ] *= k[d][ i[d] ]
        end
    end

    fourier_density = FFTW.rfft( chargedensity, box.space_dims )
    efield = Array{Array{Float64}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        efield[d] = FFTW.irfft( integrate[d] .* fourier_density, box.Nx[1], box.space_dims )
    end

    return efield
end

"""
```julia
get_electrostatic_energy( box::Box, chargedensity::Array{Float64} )
```

Obtain the electrostatic energy from the charge density.

# Notes
* The electrostatic energy is calculated in Fourier space as

  ``P_k = \\int \\rho_k \\Phi_k dk``,

  where ``P_k`` is the transformed electrostatic energy, ``\\rho_k`` is the
  transformed charge density, ``\\Phi_k = \\rho_k / |k|^2`` is the transformed
  electrostatic potential, and ``k`` is the Fourier-conjugate variable of ``x``.

* The `chargedensity` may depend upon time on its last axis, in which case
  the result will also depend upon time in the same manner.

"""
function get_electrostatic_energy( box::Box, chargedensity::Array{Float64} )

    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0

    es = abs2.( FFTW.rfft( chargedensity, box.space_dims ) ) ./ k2

    N = size( es )

    # Need to double x-wavenumbers since rfft is used
    rescale_axis = [ i == 1 ? 1 : 1:N[i] for i in 1:length(N) ]
    es[rescale_axis...] *= 0.5

    es = reducedims(es, dims = box.space_dims )

    return (prod(box.dx) / prod(box.Nx) ) * es
end

"""
```julia
get_power_per_mode( box::Box, chargedensity::Array{Float64} )
```

Obtain the power spectrum of the electrostatic energy in space.

# Notes
* If the `chargedensity` depends on time, the result will also depend on time.
"""
function get_power_per_mode( box::Box, chargedensity::Array{Float64} )

    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0

    es = abs2.( FFTW.rfft( chargedensity, box.space_dims )) ./ k2
    return (prod(box.dx) / prod(box.Nx)) * es # TODO: check normalizations
end



"""
```julia
get_dispersion_relation(box::Box, chargedensity::Array{Float64})
```

Obtain the energy density Fourier-transformed in space and time.

# Notes
* The transformation in time is performed with a `FFTW.fft`,
  and the frequency may be obtained through [`wavevector`](@ref Vlasova.wavevector).
"""
function get_dispersion_relation(box::Box, chargedensity::Array{Float64})

    k2 = get_k2(box); k2[1] = Inf # So that the inverse yields 0.0

    disprel = abs2.( FFTW.rfft( chargedensity ) ) ./ k2

    return disprel
end

# function get_electric_field(;
#                             box::Box, potential::Array{Float64})

#     Nx2p1, fourier_axis = get_rfft_dims( box.x )
#     k = rfft_wavevector( box.x )

#     integrate = Array{Array{Complex{Float64}}}(undef, box.number_of_dims)
#     for d in 1:box.number_of_dims
#         integrate[d] = ones(fourier_axis)
#         for i in fourier_axis
#             integrate[d][ i ] *= k[d][ i[d] ]
#         end
#     end

#     phik = FFTW.rfft( potential, box.space_dims )
#     efield = Array{Array{Float64}}(undef, box.number_of_dims)
#     for d in 1:box.number_of_dims
#         efield[d] = FFTW.irfft( integrate[d] .* phik, box.Nx[1], box.space_dims )
#     end

#     return efield
# end


# """
#     Generate an electrostatic potential with the size of the Box `box` at
#     an instant `time`.

#     The potential is the superposition of 1d potentials of the form
#         ``\\Phi = \\sum_i \\Phi_i \\cos ( k_i x_i - \\omega_i t)``
#     where ``i = x, y, z``...
# """
# function get_potential(box::Box, time::Real; amplitude, wavevector, frequency, time_integrated = false)
#     # TODO: This is slow
#     same_length = box.number_of_dims == length(amplitude) == length(wavevector) == length(frequency)
#     @assert same_length "The amplitude, wavevector and frequency must be of the size of the box provided"

#     if time_integrated
#         pot =  [ sum( -amplitude[d] * sin.( wavevector[d] * box.x[d][i[d]]
#                                             - frequency[d] * time  ) / frequency[d]
#                       for d in box.dim_axis )
#                  for i in CartesianIndices( box.Nx ) ]
#     else
#         pot =  [ sum( amplitude[d] * cos.( wavevector[d] * box.x[d][i[d]]
#                                            - frequency[d] * time  )
#                       for d in box.dim_axis )
#                  for i in CartesianIndices( box.Nx ) ]
#     end

#     return pot
# end
