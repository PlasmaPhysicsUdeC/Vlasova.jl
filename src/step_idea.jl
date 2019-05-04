# Integrator coefficients

function integrator()
    
    position = [ 1/3, 1/3, 1/3 ]
    velocity = [ 1/2, 1/2 ]

    position_first = true
    fuse_advections = true

    return position, velocity, position_first, fuse_advections
end

function integration_coefficients( dt )
    pos, vel, is_x_first, fuse = verlet()
    is_x_shorter = length(pos) < length(vel)
    is_v_shorter = length(pos) > length(vel)
    is_valid = (is_x_first && !is_x_shorter) || (!is_x_first && !is_v_shorter)
    pos_len, vel_len = length(pos), length(vel)
    @assert (abs( pos_len - vel_len) ) <= 1) "There are $pos_len coefficients for the space advection, and $vel_len for the velocity advections"
    @assert is_valid "In the integrator, the first advection can't have less coeffs. than the second."
    if fuse
        if is_x_first
            position[end] += position[1]
        else
            velocity[end] += velocity[1]
        end
    end
    position *= dt
    velocity *= dt

    return position, velocity
end

function step!(plasma, chargedensity, electricfield,
               poisson!, space_advection!, velocity_advection!;
               stepnumber::Integer, velocity_filtering::Bool)

    pos, vel, position_first, fuse = integrator()
    
    if position_first
        if fuse
            for i in 1:(length(pos)-1)
                velocity_advection!( plasma, electricfield, advection_number = i)
                
                space_advection!( plasma, advection_number = i + 1)
                
                # Obtain electric field
                get_density!( chargedensity, plasma )
                poisson!( electricfield, chargedensity )
            end
        else
            for i in 1:length(pos)
                space_advection!( plasma, advection_number = i)
                
                # Obtain electric field
                get_density!( chargedensity, plasma )
                poisson!( electricfield, chargedensity )

                velocity_advection!( plasma, electricfield, advection_number = i)
            end
        end
    else
        if fuse
            for i in 1:(length(vel)-1)
                space_advection!( plasma, advection_number = i + 1)
                
                # Obtain electric field
                get_density!( chargedensity, plasma )
                poisson!( electricfield, chargedensity )

                velocity_advection!( plasma, electricfield, advection_number = i)
            end
        else
            for i in 1:length(vel)
                velocity_advection!( plasma, electricfield, advection_number = i)
                
                space_advection!( plasma, advection_number = i)
                
                # Obtain electric field
                get_density!( chargedensity, plasma )
                poisson!( electricfield, chargedensity )
            end
        end
    end
    return 0;
end
