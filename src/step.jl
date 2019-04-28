# Integrator coefficients
function integration_coefficients( dt )
    position = -1im * [ dt/2 ]
    velocity =  1im * [  dt  ]

    return position, velocity
end

function step!(plasma, chargedensity, electricfield,
               poisson!, space_advection!, velocity_advection!;
               stepnumber::Integer, velocity_filtering::Bool)
    # Integrator BAB from Omelyan (2003)

    # First space advection
    space_advection!( plasma, advection_number = 1)
        
    # Obtain electric field
    get_density!( chargedensity, plasma )
    poisson!( electricfield, chargedensity )

    # Last velocity advection (+ First one)
    
    velocity_advection!( plasma, electricfield,
                         advection_number = 1,
                         advect_twice = true,
                         filtering = velocity_filtering )
    return 0;
end
