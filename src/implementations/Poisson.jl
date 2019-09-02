@doc raw"""
    Solve the Poisson equation to find the electric field generated
    by a charge distribution:

    ```math
       \vec E = - \nabla \Phi
       \nabla^2 \Phi = \rho
    ```
    where ``\vec E`` is the electric field and ``\rho`` is the charge density.

    When called with the electricfield and the chargedensity it works in place, but
    if its called with just the chargedensity, it returns the electricfield.
"""
function (poisson::Poisson)(electricfield, chargedensity; external_potential = Zero() )
    LinearAlgebra.mul!( poisson.fourier_density, poisson.plan, chargedensity)

    # Include external potential #TODO: fix this
    if typeof(external_potential) != Zero
       poisson.fourier_density .+= (poisson.plan * external_potential) .* poisson.pot2dens
    end

    for d in axes(electricfield, 1)
        @views LinearAlgebra.ldiv!( electricfield[d], poisson.plan, poisson.dens2field[d] .* poisson.fourier_density )
    end

    return nothing;
end

"""
    Get the electricfield from chargedensity
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
    Obtain the gradient forcing therm `` \vec g = \nabla |\vec E|^2 `` to use gradient integrators
"""
function get_gradient_correction!(grad, p::Poisson, E)
    Ndims = length(E)
    @views @. grad[1] = abs2( E[1] )
    for i in 2:Ndims
        @views @. grad[1] += abs2( E[i] )
    end

    @views LinearAlgebra.mul!(p.fourier_density, p.plan, grad[1])

    for i in 1:Ndims
        @views LinearAlgebra.ldiv!(grad[i], p.plan, -1im .* p.k[i] .* p.fourier_density)
    end

    return nothing;
end
