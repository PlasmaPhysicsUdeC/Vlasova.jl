"""
    Solve the Poisson equation to find the electric field generated
    by a charge distribution:
    
       \\nabla \\cdot \\vec E = -4\\pi \\rho
    where \\vec E is the electric field and \\rho is the charge density.

    When called with the electricfield and the chargedensity it works in place, but
    if its called with just the chargedensity, it returns the electricfield.
"""
function (poisson::Poisson)(electricfield, chargedensity; external_potential = Zero() )
    LinearAlgebra.mul!( poisson.fourier_density, poisson.plan, chargedensity)
    
    if typeof(external_potential) != Zero
       poisson.fourier_density .+= (poisson.plan * external_potential) .* poisson.pot2dens
    end
    
    for d in axes(electricfield, 1)
        LinearAlgebra.ldiv!( electricfield[d], poisson.plan, poisson.dens2field[d] .* poisson.fourier_density )
    end
    return 0;
end

"""
    Takes only the chargedensity and returns the electric field
"""
function (poisson::Poisson)(chargedensity; external_potential = Zero())
    Nx = size( chargedensity )
    electricfield = Array{Array{Float64}}(undef, length(Nx))
    for d in axes(electricfield, 1)
        electricfield[d] = similar(chargedensity)
    end
    poisson( electricfield, chargedensity, external_potential = external_potential) # Call in-place version
    return electricfield;
end
