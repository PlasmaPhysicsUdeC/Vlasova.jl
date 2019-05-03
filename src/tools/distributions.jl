export maxwellian1d, maxwellian2d, twostream1d, bump_on_tail1d

# 1D normalized Maxwellian
function maxwellian1d(v; vt::Float64 = 1.0, vd = 0.0::Float64) # Velocity may be anything (Float64, Array{Float64}, etc...)
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: vt = 1.0, vd = 0.0
    =#

    c = 1.0/(vt*sqrt(2pi)) # Normalization constant
    return @. c * exp( -0.5*( (v - vd)/vt )^2 )
end

function maxwellian2d(vx, vy; vtx::Float64 = 1.0, vty::Float64 = 1.0, vdx::Float64 = 0.0, vdy::Float64 = 0.0)
    #=
    Returns a normalized 2-dimensional Maxwellian.
    If not provided: Vtx, Vty = 1.0 and v0x, v0y = 0.0
    =#
    Mx = maxwellian1d(vx, vt = vtx, vd = vdx)
    My = maxwellian1d(vy, vt = vty, vd = vdy)
    
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
    c1 = nb / (vtb * sqrt(2pi)) # Normalization constant
    return @. c1 * exp( -0.5*v^2 ) + c2 * exp( -0.5 * ((v - vdb)/vtb)^2 )
end
