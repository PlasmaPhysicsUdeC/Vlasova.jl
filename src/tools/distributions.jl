export maxwellian1d,
    maxwellian2d,
    twostream1d,
    bump_on_tail1d,
    bgk1d

"""
```julia
maxwellian1d(v; vt::Real = 1.0, vd = 0.0::Real)
```

Generate a 1-d normalized Maxwellian with thermal velocity `vt` and drift velocity `vd`.
"""
function maxwellian1d(v; vt::Real = 1.0, vd = 0.0::Real)

    c = 1.0/(vt*sqrt(2pi)) # Normalization constant
    return @. c * exp( -0.5*( (v - vd)/vt )^2 )
end

"""
```julia
maxwellian1d(box::Box; vt::Real = 1.0, vd = 0.0::Real, dim::Int = 1)
```

Generate a 1-d normalized Maxwellian from a [`Vlasova.Box`](@ref) element.

# Notes
* Since the distribution is 1-dimensional, the keyword `dim` is available to use
the velocity along the dimension `dim` of the `box`.
"""
@inline function maxwellian1d(box::Box; vt::Real = 1.0, vd = 0.0::Real, dim::Int = 1)

    return maxwellian1d(box.v[dim], vt = vt, vd = vd)
end

"""
```julia
maxwellian2d(vx, vy; vtx::Real = 1.0, vty::Real = 1.0, vdx::Real = 0.0, vdy::Real = 0.0)
```

Generate a 2-d normalized Maxwellian.

# Notes
* `vtx` and `vty` are the thermal velocities along the first and second dimension, respectively.
* `vdx` and `vdy` are the drift velocities along the first and second dimension, respectively.
"""
function maxwellian2d(vx, vy; vtx::Real = 1.0, vty::Real = 1.0, vdx::Real = 0.0, vdy::Real = 0.0)
    Mx = maxwellian1d(vx, vt = vtx, vd = vdx)
    My = maxwellian1d(vy, vt = vty, vd = vdy)

    return Mx*My'
end

"""
```julia
maxwellian2d(box::Box; vtx::Real = 1.0, vty::Real = 1.0, vdx::Real = 0.0, vdy::Real = 0.0)
```

Generate a 2-d normalized Maxwellian from a [`Vlasova.Box`](@ref) element.

# Notes
* This function is intended to work with a 2-dimensional box. In other case, it will trow error.
* `vtx` and `vty` are the thermal velocities along the first and second dimension, respectively.
* `vdx` and `vdy` are the drift velocities along the first and second dimension, respectively.
"""
function maxwellian2d(box::Box; vtx::Real = 1.0, vty::Real = 1.0, vdx::Real = 0.0, vdy::Real = 0.0)
    @assert (box.number_of_dims == 2) "The box provided to this function must be 2-dimensional."

    Mx = maxwellian1d(box.v[1], vt = vtx, vd = vdx)
    My = maxwellian1d(box.v[2], vt = vty, vd = vdy)

    return Mx*My'
end

"""
```julia
twostream1d(v; vt1::Real = 1.0, vt2::Real = 1.0, vd1::Real = -2.0, vd2::Real = 2.0, n1::Real = 0.5, n2::Real = 0.5)
```

Generate a 1-d monospecies Two Stream distribution.

# Notes
* `vt1` and `vt2` are the thermal velocities of the first and second stream, respectively.
* `vd1` and `vd2` are the drift velocities of the first and second stream, respectively.
* `n1` and `n2` are the densities of the first and second stream, respectively.
"""
function twostream1d(v; vt1::Real = 1.0, vt2::Real = 1.0, vd1::Real = -2.0, vd2::Real = 2.0, n1::Real = 0.5, n2::Real = 0.5)
    @assert isapprox(n1 + n2, 1.0) "The densities of the two streams do not add up to 1"
    # Normalization constant
    c1 = n1 / ( vt1*sqrt(2pi) )
    c2 = n2 / ( vt2*sqrt(2pi) )
    return @. c1 * exp( -0.5*( (v - vd1)/vt1 )^2 ) + c2 * exp( -0.5*( (v - vd2)/vt2 )^2 )
end

"""
```julia
bump_on_tail1d(v; vtb::Real = 0.5, vdb::Real = 4, nc::Real = 0.9, nb::Real = 0.1)
```
Generate a 1-d Bump on Tail distribution.

# Notes
* The thermal velocity of the core is 1.
* `vtb` is the thermal velocity of the beam.
* `vdb` is the drift velocity of the beam.
* `nc` and `nb` are the densities of the core and the beam, respectively.
"""
function bump_on_tail1d(v; vtb::Real = 0.5, vdb::Real = 4, nc::Real = 0.9, nb::Real = 0.1)
    @assert isapprox(nc + nb, 1.0) "The densities of the core and the beam do not add up to 1"
    # Normalization constants
    c1 = nc / sqrt(2pi)
    c2 = nb / (vtb * sqrt(2pi))
    return @. c1 * exp( -0.5 * v^2 ) + c2 * exp( -0.5 * ((v - vdb)/vtb)^2 )
end


"""
```julia
bgk1d( box; amplitude::Real, wavenumber::Real, vphi::Real, dim::Int = 1 )
```

Generate a 1-dimensional BGK state.

# Notes

* The BGK state is constructed according to the mathematical description given by [`Silantyev et al. (2017-1)`](https://aip.scitation.org/doi/10.1063/1.4979289), derived from the potential ``\\Phi(x) = -\\phi_0 \\cos( k x ),`` and a Maxwellian equilibrium distribution.

* If a 2-d box is provided, the keyword `dim` may be used to select the dimension of the box used to construct the BGK state.
"""
function bgk1d( box; amplitude::Real, wavenumber::Real, vphi::Real, dim::Int = 1 ) # Use FastGaussQuadrature instead.
    @assert box.Lx[dim] == 2pi / wavenumber "The wavenumber provided does not fit exactly once in the box"

    # Consider the 1-d phase space along dimension dim.
    x = box.x[dim]
    v = box.v[dim]

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
                a = acos( -W / amplitude ) / wavenumber - 1e-12 # TODO: Try to avoid this later.

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
