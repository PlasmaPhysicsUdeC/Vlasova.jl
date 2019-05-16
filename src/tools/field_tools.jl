export get_kinetic_energy,
    get_density,
    get_density!,
    get_k2,
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
    Integrates a distribution function in the velocity space to obtain the charge density

    Requires:
    * distribution: Array of Float64
    * box: Element of type Box
    
    Optional, keyword:
    * set_zero_mean: [false] Substract the mean of the resulting chargedensity to ensure quasineutrality
    
    Returns:
    * chargedensity: Array of Float64
"""
function get_density(distribution::Array{Float64}, box::Box;
                     set_zero_mean::Bool = false)
    # Integrate
    chargedensity = -prod(box.dv) * dropdims( sum(distribution, dims = box.velocity_dims ), dims = box.velocity_dims )
    # Quasineutrality
    set_zero_mean ? (chargedensity .= chargedensity .- mean(chargedensity )) : nothing
    return chargedensity
end

"""
    Integrates a distribution function in the velocity space to obtain the charge density in place!
    ** This version in no more efficient than it's not-in-place counterpart! **

    Requires:
    * chargedensity: Array of Float64
    * distribution: Array of Float64
    * box: Element of type Box
    
    Optional, keyword:
    * set_zero_mean: [false] Substract the mean of the resulting chargedensity to ensure quasineutrality
    
    Returns:
    * Error code 0
"""
function get_density!(chargedensity::Array{Float64}, distribution::Array{Float64}, box::Box;
                      set_zero_mean::Bool = false)
    
    chargedensity .= get_density( distribution, box, set_zero_mean = set_zero_mean)
    return 0;
end

"""
    Obtain the squared wavevector,

    `k^2 = k_1^2 + k_2^2 + ... k_n^2`

    where k_1 is a half of the first wavevector as it would be used on a real DFT.
"""
function get_k2( box::Box )
    Nx2p1, fourier_axis = get_rfft_dims( box )

    k = rfft_wavevector( box.x )

    k2 = zeros( Nx2p1 )
    for d in 1:box.number_of_dims, i in fourier_axis
        k2[i] += ( k[d][ i[d] ] )^2
    end

    return k2
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
function get_electric_field(chargedensity::Array{Float64}, box::Box) # TODO: check!
    
    Nx2p1, fourier_axis = get_rfft_dims( box )
    
    k = rfft_wavevector( box.x )
    k2 = get_k2( box ); k2[1] = Inf  # So that the inverse yields 0.0
    
    integrate = Array{Array{Complex{Float64}}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        integrate[d] = -1im ./ k2
        for i in fourier_axis
            integrate[d][ i ] *= k[d][ i[d] ]
        end
    end
    
    fourier_density = FFTW.rfft( chargedensity )
    efield = Array{Array{Float64}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        efield[d] = FFTW.irfft( integrate[d] .* fourier_density, box.Nx[1], box.space_dims )
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
