
function (vadv::VelocityAdvection)(plasma::Plasma, electricfield, grad, prop;
                                   advection_number::Int64,
                                   gradient_number::Int64,
                                   is_gradient_advection::Bool,
                                   filtering::Bool)
    
    for s in 1:plasma.number_of_species
        
        # DF to velocity Fourier-space
        LinearAlgebra.mul!(vadv.transformed_DF, vadv.plan, plasma.species[s].distribution )
        
        # Apply high-freq filter
        filtering ? lowpass_velocityfilter!( vadv.transformed_DF, vadv.filter, plasma.box.number_of_dims ) : nothing

        # Prepare reduced propagator
        # Electricfield coefficient
        coef = vadv.advection_coefficients[ advection_number ][ s ]
        
        dp = [vadv.wavevector[d][2] - vadv.wavevector[d][1] # delta v_wavevector
              for d in 1:plasma.box.number_of_dims ]
        
        if is_gradient_advection
            grad_coef = vadv.gradient_coefficients[ gradient_number ][ s ]
            
            # Propagator dependency on E with gradient correction # TODO: Not working yet
            for d in 1:plasma.box.number_of_dims
                @. prop[d] = cis( coef * ( electricfield[d] + grad_coef * grad[d] ) * dp[d] )
            end
        else
            # Propagator dependency on E
            for d in 1:plasma.box.number_of_dims
                @. prop[d] = cis( coef * electricfield[d] * dp[d] )
            end            
        end

        # Advect with corrected forcing
        _velocity_advection!( vadv.transformed_DF, prop, size(vadv.transformed_DF)...)

        # Return to real space
        LinearAlgebra.ldiv!( plasma.species[s].distribution, vadv.plan, vadv.transformed_DF )
        
    end
    return 0;
end

function lowpass_velocityfilter!(transformed_DF, filter, dims)
    
    if dims == 1
        for i in CartesianIndices(filter)
            @views @. transformed_DF[:, i] *= filter[i]
        end
    else
        for i in CartesianIndices(filter)
            @views @. transformed_DF[:, :, i] *= filter[i]
        end
    end
    
    return 0;
end


# 1D
function _velocity_advection!(transformed_DF, prop, Nx, Nvx)
    # Calculate propagator as exp(coef * Ex * dux)^j, where ux = dux*j and j = 0,1,2...
    e0 = ones(Complex{Float64}, Nx )
    for i in 1:Nvx
        @views @. transformed_DF[:, i] *= e0
        @views @. e0 *= prop[1]
    end
    
    return 0;
end

# 2D
function _velocity_advection!(transformed_DF, prop, Nx, Ny, Nvx, Nvy)
    
    ex = Array{Complex{Float64}}(undef, Nx, Ny, Nvx)
    e0 = ones(Complex{Float64}, Nx, Ny)
    for j in 1:Nvx
        @views @. ex[:, :, j] = e0
        @views @. e0 *= prop[1]
    end

    Nvy2 = div(Nvy, 2)    
    fill!(e0, one(Complex{Float64}))
    for i in 1:Nvy2 # p_y Index growing
        for j in 1:Nvx
            @views @. transformed_DF[:, :, j, i] *= e0 * ex[:, :, j]
        end
        @views @. e0 *= prop[2]
    end

    @. prop[2] = 1 / prop[2]    # momentum indices are now negative
    @. e0 =  prop[2]
    for i in Nvy:-1:(Nvy2+1) # p_y Index decreasing
        for j in 1:Nvx
            @views @. transformed_DF[:, :, j, i] *= e0 * ex[:, :, j]
        end
        @views @. e0 *= prop[2]
    end
    
    return 0;
end
