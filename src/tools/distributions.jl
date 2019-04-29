export maxwellian1d, maxwellian2d

# 1D normalized Maxwellian
function maxwellian1d(velocity; thermal_velocity::Float64 = 1.0, drift_velocity = 0.0::Float64) # Velocity may be anything (Float64, Array{Float64}, etc...)
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: thermal_velocity = 1.0, drift_velocity = 0.0
    =#

    normalization = 1.0/(thermal_velocity*sqrt(2pi)) # Normalization constant
    return @. normalization * exp( -0.5*( (velocity - drift_velocity)/thermal_velocity )^2 )
end

function maxwellian2d(vx, vy; vtx::Float64 = 1.0, vty::Float64 = 1.0, v0x::Float64 = 0.0, v0y::Float64 = 0.0)::Array{Float64, 2}
    #=
    Returns a normalized 2-dimensional Maxwellian.
    If not provided: Vtx, Vty = 1.0 and v0x, v0y = 0.0
    =#
    Mx = maxwellian1d(vx, thermal_velocity = vtx, drift_velocity = v0x)
    My = maxwellian1d(vy, thermal_velocity = vty, drift_velocity = v0y)
    
    return Mx*My'
end
