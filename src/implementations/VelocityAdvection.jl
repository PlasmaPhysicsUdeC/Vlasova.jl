# TODO: THREADS
"""
Apply a velocity advection over the plasma.
"""
function (vadv::VelocityAdvection)(plasma::Plasma, electricfield, grad, prop;
                                   advection_number::Int64,
                                   gradient_number::Int64,
                                   is_gradient_advection::Bool,
                                   filtering::Bool)

    for s in 1:plasma.number_of_species

        # DF to velocity Fourier-space
        TimerOutputs.@timeit_debug timer "Fourier transform" begin
            LinearAlgebra.mul!(vadv.transformed_DF, vadv.plan, plasma.species[s].distribution )
        end

        # Apply high-freq filter
        TimerOutputs.@timeit_debug timer "filtering" begin
            filtering ? lowpass_velocityfilter!( vadv.transformed_DF, vadv.filter, plasma.box ) : nothing
        end

        # Prepare reduced propagator
        # Electricfield coefficient
        coef = vadv.advection_coefficients[ :, s, advection_number ]

        TimerOutputs.@timeit_debug timer "prepare reduced propagator" begin
            if is_gradient_advection
                grad_coef = vadv.gradient_coefficients[ :, s, gradient_number]

                # Propagator dependency on E with gradient correction # TODO: Not working yet
                for d in 1:plasma.box.number_of_dims
                    Strided.@strided prop[d] .= cis.( - coef[d] .* electricfield[d] .+ grad_coef[d] .* grad[d] )
                end
            else
                # Propagator dependency on E
                for d in 1:plasma.box.number_of_dims
                    Strided.@strided prop[d] .= cis.( - coef[d] .* electricfield[d] )
                end
            end
        end

        # Perform advection with corrected force
        TimerOutputs.@timeit_debug timer "calculate and apply propagator" begin
            _velocity_advection!( vadv.transformed_DF, prop, size(vadv.transformed_DF)...)
        end

        # Return to real space
        TimerOutputs.@timeit_debug timer "inverse Fourier transform" begin
            LinearAlgebra.ldiv!( plasma.species[s].distribution, vadv.plan, vadv.transformed_DF )
        end

    end
    return nothing;
end

"""
Apply element-wise product of the distribution function in the velocity-transformed space against a given filter.
"""
function lowpass_velocityfilter!(transformed_DF, filter, box)

    Threads.@threads for i in CartesianIndices(filter)
        @inbounds @views @. transformed_DF[box.space_axes..., i] *= filter[i]
    end

    return nothing;
end


"""
[1D] Perform the product of the distribution in the velocity-transformed space with the corresponding exponential propagator.

The product is performed considering that

`` \\exp\\left( -i  E_x[j]  p_x[k] \\left) = \\exp \\left( - i \\Delta p_x E_x[j]  \\right)^k`` where ``k = 0, 1, 2, \\ldots, N_{vx}/2``.

"""
function _velocity_advection!(transformed_DF, prop, Nx, Nvx)
    # prop[i] = exp( -i * dp * E[i] )
    e0 = ones(Complex{Float64}, Nx ) # prop ^ j, with j = 0
    for i in 1:Nvx
        @inbounds @views @. transformed_DF[:, i] *= e0
        @inbounds @views @. e0 *= prop[1] # j = 1, 2, 3 ...
    end

    return nothing;
end

"""
[2D] Perform the product of the distribution in the velocity-transformed space with the corresponding exponential propagator.

The product is performed considering that

`` \\exp\\left[ -i  ( E_x[j]  p_x[k] + E_y[l] p_y[m] \\left] = \\exp \\left[ - i \\Delta p_x E_x[j]  \\right]^k \\exp \\left[ - i \\Delta p_y E_y[l]  \\right]^m`

 where ``k = 0, 1, 2, \\ldots, N_{vx}/2`` and  `` m = -N_{vy}/2, \\ldots, N_{vy}/2``.

"""
function _velocity_advection!(transformed_DF, prop, Nx, Ny, Nvx, Nvy)
    # Calculate the exponential propagator exp( -i * p_x * E_x ) as in the 1-d version.
    ex = Array{Complex{Float64}}(undef, Nx, Ny, Nvx)
    e0 = ones(Complex{Float64}, Nx, Ny)
    for j in 1:Nvx
        @inbounds @views @. ex[:, :, j] = e0
        @inbounds @views @. e0 *= prop[1]
    end

    Nvy2 = div(Nvy, 2)
    fill!(e0, one(Complex{Float64}))
    for i in 1:Nvy2 # p_y Index growing
        for j in 1:Nvx
            @inbounds @views @. transformed_DF[:, :, j, i] *= e0 * ex[:, :, j]
        end
        @inbounds @views @. e0 *= prop[2]
    end

    @. prop[2] = 1 / prop[2]    # momentum indices are now negative
    @. e0 =  prop[2]
    for i in Nvy:-1:(Nvy2+1) # p_y Index decreasing
        for j in 1:Nvx
            @inbounds @views @. transformed_DF[:, :, j, i] *= e0 * ex[:, :, j]
        end
        @inbounds @views @. e0 *= prop[2]
    end

    return nothing;
end
