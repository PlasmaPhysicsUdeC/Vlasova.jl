export get_electrostatic_energy, get_kinetic_energy, get_density

function get_electrostatic_energy( chargedensity::Array{Float64}, x )
    #=
    Obtaing the electrostatic energy from the charge density.
    If the charge density dependes upon time, the result will also depend upon time.

    The electrostatic energy is calculated in Fourier space as
    Energy = \int \rho_k* \Phi_k dk,
    where \rho_k is the transformed charge density, \Phi_k = rho_k / |k|^2 is the
    transformed electrostatic potential, and k is the Fourier-conjugate variable of x.
    =#
    
    Nx = size( chargedensity )
    number_of_dims = length( Nx )
    dx = Tuple( x[d][2] - x[d][1] for d in 1:number_of_dims )
    
    Nx2p1 = Tuple( (i == 1) ? fld( Nx[1], 2 ) + 1 : Nx[i] for i in 1:number_of_dims )
    fourier_axis = CartesianIndices( Nx2p1 )

    k = Array{Array{Float64, 1}}(undef, number_of_dims)
    k[1] = rfft_wavevector( x[1] )
    for d in 2:number_of_dims
        k[d] = wavevector( x[d] )
    end

    k2 = zeros( fourierDims )
    for d in dim_axis, i in fourier_axis
        k2[i] += ( k[d][ i[d] ] )^2
    end
    k2[1] = Inf  # So that the inverse yields 0.0

    es = sum( abs2.( FFTW.rfft( chargedensity, space_dims ) ./ k2 ), dims = space_dims )

    return (0.5 * prod(dx) / prod(Nx) ) * dropdims( es, dims = space_dims)
end

function get_kinetic_energy(distribution::Array{Float64}, v, dx; temperature = 1.0);
    #=
    Requires:
    * dSpace, dVelocity
    * distributionFunction[ Nx..., Nv... ]

    [optional]
    * temperature of the specie with respect to electrons

    Returns:
    * kineticEnergy (scalar)
    =#
    N = size( distribution )
    number_of_dims = fld( length( N ), 2)
    space_dims = Tuple( 1:number_of_dims )
    dv = Tuple( v[d][2] - v[d][1] for d in space_dims )
    
    kinen = dropdims( sum( distribution, dims = space_dims ), dims = space_dims )

    for d in 1:number_of_dims
        @. kinen *= @views v[d]^2
        kinen = dropdims( sum(kinen, dims = 1 ), dims = 1)
    end

    return temperature * 0.5 * prod(dx) * prod(dv) * kinen[1];
end

function get_density(distribution::Array{Float64}, dv)
    #=
    Integrates a distribution function in the velocity space to obtain the charge density.

    Requires:
    * Nx, dVelocity
    * distribution[Nx..., Nv...]

    Returns:
    * chargedensity[ Nx... ]
    =#

    N = size( distribution )
    number_of_dims = fld( length( N ), 2)
    velocity_dims = Tuple( 1:number_of_dims ) .+ number_of_dims
    
    return -prod(dv) * dropdims( sum(distribution, dims = velocity_dims ), dims = velocity_dims )
end

function get_density!(chargedensity::Array{Float64}, distribution::Array{Float64}, dv, set_zero_mean::Bool = false)
    #=
    Integrates a distribution function in the velocity space to obtain the charge density in place!.

    Requires:
    * Nx, dVelocity
    * chargedensity[Nx...]
    * distribution[Nx..., Nv...]

    Keywords:
    * If set_zero_mean = true is passed, quasineutrality is ensured. This behavior is disabled by default.

    Returns:
    * Error code 0;
    =#

    chargedensity .= get_density( distribution, dv)

    if set_zero_mean
        chargedensity .-= mean(chargedensity)
    end
    return 0;
end
