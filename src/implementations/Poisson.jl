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

"""
    Obtains the gradient forcing therm ` grad = \nabla |\vec E|^2 ` to use gradient integrators
"""
function gradient_force!(grad, p::Poisson, E)
    Ndims = length(E)
    @. grad[1] = abs2( E[1] )
    for i in 2:Ndims
        @. grad[1] += abs2( E[i] )
    end
    LinearAlgebra.mul!(p.fourier_density, p.plan, grad[1])
    # Swap real and imaginary parts outside the loop
    @. p.fourier_density *= 1im
    for i in 1:Ndims
        LinearAlgebra.ldiv!(grad[i], p.plan, p.k[i] .* p.fourier_density) # TODO: This will fail in 2d
    end
    
    return 0;
end
    
