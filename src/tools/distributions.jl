
# 1D normalized Maxwellian
function maxwellian1d(velocity::Array{Float64, 1}; thermalVelocity::Float64 = 1.0, driftVelocity = 0.0::Float64) # Args as keywords
    #=
    Returns a normalized 1-dimensional Maxwellian.
    If not provided: thermalVelocity = 1.0, driftVelocity = 0.0
    =#

    normalizationConstant = 1.0/(thermalVelocity*sqrt(2pi)) # Normalization constant
    return @. normalizationConstant * exp( -0.5*( (velocity-driftVelocity)/thermalVelocity )^2 )
end


# 1D two stream instability, monospecies
function twoStream1sp1d(velocity::Array{Float64,1}; thermalVelocity1::Float64 = 1.0, thermalVelocity2::Float64 = 1.0, driftVelocity1::Float64 = -2.0, driftVelocity2::Float64 = 2.0, nb1::Float64 = 0.5, nb2::Float64 = 0.5)
    #=
    Returns a normalized 1-dimensional two stream distribution.
    If not provided: thermalVelocity1 = thermalVelocity = 1.0, driftVelocity1 = driftVelocity2 = -2.0 
    =#

    normalizationConstant1 = nb1/(thermalVelocity1*sqrt(2pi)) # Normalization constant
    normalizationConstant2 = nb2/(thermalVelocity2*sqrt(2pi)) # Normalization constant
    return @. normalizationConstant1 * exp( -0.5*(velocity - driftVelocity1)^2/thermalVelocity1^2 ) + normalizationConstant2 * exp( -0.5*(velocity - driftVelocity2)^2/thermalVelocity2^2 )
end


# 1D bump on tail distribution, monospecies
function bumpOnTail1d(velocity::Array{Float64,1}; thermalVelocityBeam::Float64 = 0.5, driftVelocityBeam::Float64 = 4.5, densityCore::Float64 = 0.9, densityBeam::Float64 = 0.1)
    #=
    Returns a normalized 1-dimensional two stream distribution.
    If not provided: thermalVelocity1 = thermalVelocity = 1.0, driftVelocity1 = driftVelocity2 = -4.0 
    =#

    normalizationConstant1 = densityCore/sqrt(2pi) # Normalization constant
    normalizationConstant2 = densityBeam/(thermalVelocityBeam*sqrt(2pi)) # Normalization constant
    return @. normalizationConstant1 * exp( -0.5*( velocity )^2) + normalizationConstant2 * exp( -0.5*(velocity - driftVelocityBeam)^2/thermalVelocityBeam^2)
end
