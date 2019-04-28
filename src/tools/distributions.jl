export maxwellian1d, maxwellian2d

# 1D normalized Maxwellian
function maxwellian1d(velocity; thermalVelocity::Float64 = 1.0, driftVelocity = 0.0::Float64) # Velocity may be anything (Float64, Array{Float64}, etc...)
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: thermalVelocity = 1.0, driftVelocity = 0.0
    =#

    normalizationConstant = 1.0/(thermalVelocity*sqrt(2pi)) # Normalization constant
    return @. normalizationConstant * exp( -0.5*( (velocity-driftVelocity)/thermalVelocity )^2 )
end

function maxwellian2d(vx, vy; vtx::Float64 = 1.0, vty::Float64 = 1.0, v0x::Float64 = 0.0, v0y::Float64 = 0.0)::Array{Float64, 2}
    #=
    Returns a normalized 2-dimensional Maxwellian.
    If not provided: Vtx, Vty = 1.0 and v0x, v0y = 0.0
    =#
    Mx = maxwellian1d(vx, thermalVelocity = vtx, driftVelocity = v0x)
    My = maxwellian1d(vy, thermalVelocity = vty, driftVelocity = v0y)
    
    return Mx*My'
end
