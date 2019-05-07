"""
    This function solves the Poisson equation to find the electric field generated
    by a charge distribution:
    
       \\nabla \\cdot \\vec E = -4\\pi \\rho
    where \\vec E is the electric field and \\rho is the charge density.

    When called with the electricfield and the chargedensity it works in place, but
    if its called with just the chargedensity, it returns the electricfield.
"""
function (poisson::Poisson)(electricfield, chargedensity )
    LinearAlgebra.mul!( poisson.fourier_density, poisson.plan, chargedensity)
    
    for d in length(electricfield)
        LinearAlgebra.ldiv!( electricfield[d], poisson.plan, poisson.integrate[d] .* poisson.fourier_density )
    end
    return 0;
end

"""
    Takes only the chargedensity and returns the electric field
"""
function (poisson::Poisson)(chargedensity)
    Nx = size( chargedensity )
    electricfield = Array{Array{Float64}}(undef, length(Nx))
    for d in 1:length(Nx)
        electricfield[d] = similar(chargedensity)
    end
    poisson( electricfield, chargedensity) # Call in-place version
    return electricfield;
end
