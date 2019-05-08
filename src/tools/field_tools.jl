export get_electrostatic_energy, get_electric_field, get_kinetic_energy, get_density

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
    
    Nx2p1 = Tuple( (i == 1) ? fld( box.Nx[1], 2 ) + 1 : box.Nx[i] for i in 1:box.number_of_dims )
    fourier_axis = CartesianIndices( Nx2p1 )

    k = Array{Array{Float64, 1}}(undef, box.number_of_dims)
    k[1] = rfft_wavevector( box.x[1] )
    for d in 2:box.number_of_dims
        k[d] = wavevector( box.x[d] )
    end

    k2 = zeros( Nx2p1 )
    for d in 1:box.number_of_dims, i in CartesianIndices( Nx2p1 )
        k2[i] += ( k[d][ i[d] ] )^2
    end
    k2[1] = Inf  # So that the inverse yields 0.0

    es = sum( abs2.( FFTW.rfft( chargedensity, box.space_dims )) ./ k2, dims = box.space_dims )

    return (prod(box.dx) / prod(box.Nx) ) * dropdims( es, dims = box.space_dims)
end

"""
    Obtains the electric field from a charge density.

    In general, te electric field returned is an array of arrays, where the first array
    correspond to the electric field along the first dimension and so on.

    In the case where the charge density depends on time, the electric field will continue
    to be an array of arrays, where each of them will be the electric field along one dimensions
    having the same dependence on time as the charge density.
"""
function get_electric_field(chargedensity::Array{Float64}, box::Box) # TODO: check!

    Nx2p1 = Tuple( i == 1 ? fld(box.Nx[i], 2)+1 : box.Nx[i]
                   for i in 1:length(box.Nx) )
    
    fourier_axis = CartesianIndices( Nx2p1 )
    
    fourier_density = FFTW.rfft(chargedensity, box.space_dims )
    
    k = Array{Array{Float64, 1}}(undef, box.number_of_dims)

    k[1] = rfft_wavevector( box.x[1] )
    for d in 2:box.number_of_dims
        k[d] = wavevector( box.x[d] )
    end
    
    k2 = zeros( Nx2p1 )
    for d in 1:box.number_of_dims, i in fourier_axis
        k2[i] += ( k[d][ i[d] ] )^2
    end
    k2[1] = Inf  # So that the inverse yields 0.0. Ensure quasineutrality
    
    integrate = Array{Array{Complex{Float64}}}(undef, box.number_of_dims)
    integrate[1] = -1im ./ k2
    for d in 2:box.number_of_dims
        integrate[d] = integrate[1]
    end
    
    for d in box.dim_axis, i in fourier_axis
        integrate[d][ i ] *= k[d][ i[d] ]
    end

    efield = Array{Array{Float64}}(undef, box.number_of_dims)
    for d in 1:box.number_of_dims
        efield[d] = FFTW.irfft( integrate[d] .* fourier_density, box.Nx[1], box.space_dims )
    end
    
    return efield
end

"""
    Obtains the kinetic energy from a distribution function
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
