export maxwellian1d, maxwellian2d, twostream1d, bump_on_tail1d, bgk1d

"""
    Generate a 1-d normalized Maxwellian
    """
function maxwellian1d(v; vt::Float64 = 1.0, vd = 0.0::Float64) # Velocity may be anything (Float64, Array{Float64}, etc...)
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: vt = 1.0, vd = 0.0
    =#

    c = 1.0/(vt*sqrt(2pi)) # Normalization constant
    return @. c * exp( -0.5*( (v - vd)/vt )^2 )
end

"""
    Generate a 1-d normalized maxwellian from a Vlasova.Box element
    """
function maxwellian1d(box::Box; vt::Float64 = 1.0, vd = 0.0::Float64) # Velocity may be anything (Float64, Array{Float64}, etc...)
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: vt = 1.0, vd = 0.0
    =#

    c = 1.0/(vt*sqrt(2pi)) # Normalization constant
    return @. c * exp( -0.5*( (box.v[1] - vd)/vt )^2 )
end

"""
    Generate a 2-d normalized Maxwellian
    """
function maxwellian2d(vx, vy; vtx::Float64 = 1.0, vty::Float64 = 1.0, vdx::Float64 = 0.0, vdy::Float64 = 0.0)
    #=
    Returns a normalized 2-dimensional Maxwellian.
    If not provided: vtx, vty = 1.0 and v0x, v0y = 0.0
    =#
    Mx = maxwellian1d(vx, vt = vtx, vd = vdx)
    My = maxwellian1d(vy, vt = vty, vd = vdy)

    return Mx*My'
end

"""
    Generate a 2-d normalized Maxwellian from a Vlasova.Box element
    """
function maxwellian2d(box::Box; vtx::Float64 = 1.0, vty::Float64 = 1.0, vdx::Float64 = 0.0, vdy::Float64 = 0.0)
    #=
    Returns a normalized 2-dimensional Maxwellian.
    If not provided: vtx, vty = 1.0 and v0x, v0y = 0.0
    =#
    Mx = maxwellian1d(box.v[1], vt = vtx, vd = vdx)
    My = maxwellian1d(box.v[2], vt = vty, vd = vdy)

    return Mx*My'
end

# 1D two stream instability, monospecies
function twostream1d(v; vt1::Float64 = 1.0, vt2::Float64 = 1.0, vd1::Float64 = -2.0, vd2::Float64 = 2.0, n1::Float64 = 0.5, n2::Float64 = 0.5)
    #=
    Returns a normalized 1-dimensional two stream distribution.
    =#

    # Normalization constants
    c1 = n1 / ( vt1*sqrt(2pi) )
    c2 = n2 / ( vt2*sqrt(2pi) )
    return @. c1 * exp( -0.5*( (v - vd1)/vt1 )^2 ) + c2 * exp( -0.5*( (v - vd2)/vt2 )^2 )
end


# 1D bump on tail distribution, monospecies
function bump_on_tail1d(v; vtb::Float64 = 0.5, vdb::Float64 = 4.5, nc::Float64 = 0.9, nb::Float64 = 0.1)
    #=
    Returns a normalized 1-dimensional bump on tail distribution.
    =#

    c1 = nc / sqrt(2pi) # Normalization constant
    c2 = nb / (vtb * sqrt(2pi)) # Normalization constant
    return @. c1 * exp( -0.5*v^2 ) + c2 * exp( -0.5 * ((v - vdb)/vtb)^2 )
end


"""
    Return a 1-dimensional BGK state. Only valid for the first mode.
"""
function bgk1d( box; amplitude, wavenumber, vphi, along_dim = 1 ) # TODO: Requires QuadGK
    @assert box.Lx[along_dim] == 2pi / wavenumber "The wavenumber provided does not fit once in the box"

    x = box.x[along_dim]
    v = box.v[along_dim]

    dx = x[2] - x[1]

    Nx = size(x, 1)

    Lx = dx * Nx
    Nv = size(v, 1)

    # Functions
    vx(x, W) = sqrt( 2 * ( W + amplitude * cos( wavenumber * x) ) )
    f0(v) = exp(-0.5v^2) / sqrt(2pi)
    fpm(vx) = ( f0(vphi + vx ) + f0(vphi - vx) ) / vx
    fp(vx) = f0(vphi + vx ) / vx
    fm(vx) = f0(vphi - vx ) / vx

    # Distribution
    distribution = zeros(Nx, Nv)
    energy = similar( distribution )
    for j in 1:Nv
        for i in 1:Nx
            # Single particle energy in the wave frame
            W = 0.5 * (v[j] - vphi)^2 - amplitude * cos( wavenumber * x[i] )

            if abs( W ) < amplitude # Trapped particles
                # Integration limit + small displacement to avoid singularity
                a = acos( -W / amplitude ) / wavenumber - 1e-12
                #b = Lx - a

                # Closed path: Contribution over and under vphi
                BGKT, = QuadGK.quadgk( x -> fpm(vx(x, W)), 0, a)
                # Orbit length
                T,  = 2 .* QuadGK.quadgk(x -> 1/vx(x, W), 0, a)
            else
                # Free particles: Integration along all positions
                T, = QuadGK.quadgk( x -> 1/vx(x, W), 0, Lx )
                if v[j] > vphi # Free orbits over vphi
                    BGKT, = QuadGK.quadgk( x -> fp(vx(x, W)), 0, Lx )
                else        # Free orbits under vphi
                    BGKT, = QuadGK.quadgk( x -> fm(vx(x, W)), 0, Lx )
                end
            end
            # Path contribution over path length
            distribution[i, j] = BGKT / T
            energy[i, j] = W
        end
    end

    return distribution, energy
end
