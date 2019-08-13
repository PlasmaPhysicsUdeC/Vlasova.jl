export get_kinetic_energy,
    get_density,
    get_density!,
    get_electric_field,
    get_electrostatic_energy,
    get_power_per_mode,
    get_dispersion_relation


"""
    Obtain the kinetic energy from a distribution function
    Requires:
    * distribution: Array of Float64
    * box: Element of type Box

    Optional, keyword:
    * temperature: [1.0] The temperature normalized to electron's temperature

    Returns:
    * kinetic_energy: Float64
"""
function get_kinetic_energy(distribution::Array{Float64}, box::Box;
                            temperature::Real = 1.0);

    kinen = dropdims( sum( distribution, dims = box.space_dims ), dims = box.space_dims )

    for d in 1:box.number_of_dims
        @. kinen *= @views box.v[d]^2
        kinen = dropdims( sum( kinen, dims = 1 ), dims = 1)
    end

    return temperature * 0.5 * prod(box.dx) * prod(box.dv) * kinen;
end

"""
    Obtain the electric field from a charge distribution.

    In general, the electric field will be an array where the n-th component is the electric field
    along the n-th dimension.

    If the charge density depends upon time, the electrif field will also depend on time.

    ```@example
    using Vlasova
    box = Box("sim", 64, 64, 2pi, -6, 6 );
    chargedensity = sin.( box.x[1] );
    Ex = get_electric_field(chargedensity, box)
    ```
"""
function get_electric_field(chargedensity::Array{Float64}, box::Box)

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

function get_electric_field(;
                            potential::Array{Float64}, box::Box)

    Nx2p1, fourier_axis = get_rfft_dims( box.x )
    k = rfft_wavevector( box.x )

    integrate = Array{Array{Complex{Float64}}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        integrate[d] = ones(fourier_axis)
        for i in fourier_axis
            integrate[d][ i ] *= k[d][ i[d] ]
        end
    end

    phik = FFTW.rfft( potential, box.space_dims )
    efield = Array{Array{Float64}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        efield[d] = FFTW.irfft( integrate[d] .* phik, box.Nx[1], box.space_dims )
    end

    return efield
end

raw"""
    Obtains the electrostatic energy from the charge density.
    If the charge density dependes upon time, the result will also depend upon time.

    The electrostatic energy is calculated in Fourier space as
    Energy = \int \rho_k* \Phi_k dk,
    where \rho_k is the transformed charge density, \Phi_k = rho_k / |k|^2 is the
    transformed electrostatic potential, and k is the Fourier-conjugate variable of x.

    The required variables are:
    * chargedensity: Array of Float64
    * box: An element of type Box.
"""
function get_electrostatic_energy( chargedensity::Array{Float64}, box::Box )

    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0

    es = abs2.( FFTW.rfft( chargedensity, box.space_dims ) ) ./ k2

    N = size( es )

    rescale_axis = [ i == 1 ? 1 : 1:N[i] for i in 1:length(N) ]

    es[rescale_axis...] *= 0.5

    es = sum(es, dims = box.space_dims )

    return (prod(box.dx) / prod(box.Nx) ) * dropdims( es, dims = box.space_dims)
end

"""
Obtain the power spectrum of the electrostatic energy in space.

Depending on chargedensity, the result may depend on time.
"""
function get_power_per_mode( chargedensity::Array{Float64}, box::Box )

    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0

    es = abs2.( FFTW.rfft( chargedensity, box.space_dims )) ./ k2
    return (prod(box.dx) / prod(box.Nx)) * es # TODO: check normalizations
end



"""
    Obtain the energy density fourier-transformed in space and time.
"""
function get_dispersion_relation(chargedensity::Array{Float64}, box::Box)

    k2 = get_k2(box); k2[1] = Inf # So that the inverse yields 0.0

    disprel = abs2.( FFTW.rfft( chargedensity ) ) ./ k2

    return disprel
end
#
# TODO: This 2 functions give different results. Why?
#
function get_dispersion_relation2(chargedensity::Array{Float64}, box::Box)
    efield = get_electric_field(chargedensity, box)

    disprel = abs2.( FFTW.rfft(efield[1]) )
    for d in 2:box.number_of_dims
        disprel += abs2.( FFTW.rfft( efield[d] ) )
    end
    return disprel              # TODO: scaling?
end

"""
    Generate an electrostatic potential with the size of the Box `box` at
    an instant `time`.

    The potential is the superposition of 1d potentials of the form
     `\\Phi = \\sum_i \\Phi_i \\cos ( k_i x_i - \\omega_i t)`
    where i = x, y, z...
"""
function get_potential(box::Box, time::Real; amplitude, wavevector, frequency, time_integrated = false)
    # TODO: This is slow
    same_length = box.number_of_dims == length(amplitude) == length(wavevector) == length(frequency)
    @assert same_length "The amplitude, wavevector and frequency must be of the size of the box provided"

    if time_integrated
        pot =  [ sum( -amplitude[d] * sin.( wavevector[d] * box.x[d][i[d]]
                                            - frequency[d] * time  ) / frequency[d]
                      for d in box.dim_axis )
                 for i in CartesianIndices( box.Nx ) ]
    else
        pot =  [ sum( amplitude[d] * cos.( wavevector[d] * box.x[d][i[d]]
                                           - frequency[d] * time  )
                      for d in box.dim_axis )
                 for i in CartesianIndices( box.Nx ) ]
    end

    return pot
end
