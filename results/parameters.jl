using Vlasova

# Name to save the data
simulation_name = "huehu"

# Space nodes
Nx = (16, 16)
Nv = (32, 32)

# Space lengths
Lx = (5pi, 5pi)
vMin = (-6.0, -6.0 )            # Minimum velocities
vMax = (6.0, 6.0 )            # Maximum velocities

# Final conditions
dt = 1e-1
final_time = 100                # In electron plasma periods

# Multi-threding
FFTW_NUM_THREADS = 4
OMP_NUM_THREADS  = 4

# Multi-specie support
name = ["electrons", "protons"]
charge = [-1.0, 1.0]
mass = [1.0, 1836.15267389]
temperature = [1.0, 1.0]
perturbed = [true, false]

# Starting conditions
function perturbate!(distribution::Array{Float64}, box::Box)
    # Perturbation parameters
    mode = (1, 0)
    amplitude = 1e-3 .* (1, 0)
    
    # Perturbation of the form ( Ax*cos(kx*x) + Ay*cos(ky*y) ... )
    k_mode = @. 2pi * (mode / box.Lx) * box.x
    perturbation = zeros( box.Nx )
    for d in box.dim_axis, i in box.space_axis
        perturbation[i] += amplitude[d] * cos( k_mode[d][ i[d] ] )
    end
    
    # Ensure that perturbation preserves mass
    perturbation = perturbation .-  mean(perturbation) 

    # Apply perturbation
    for i in box.space_axis, j in box.velocity_axis
        distribution[i, j] = ( 1 + perturbation[i] ) * distribution[i, j]
    end
    return 0;
end

function initial_distribution(box::Box; perturbate::Bool = false)

    distribution = ones( box.N )

    # TODO: check this
    for d in box.dim_axis, j in box.velocity_axis, i in box.space_axis
        distribution[i, j] *= Vlasova.maxwellian1d( box.v[d][ j[d] ], driftVelocity = 0.0  )
    end

    # Apply perturbation ?
    perturbate ? perturbate!( distribution, box ) : nothing
    return distribution
end

# Create plasma box
box = Box( simulation_name, Nx, Nv, Lx, vMin, vMax )

# Create specie array
species = [ Specie(name[s],
                   charge[s],
                   mass[s],
                   temperature[s],
                   initial_distribution(box, perturbate = perturbed[s]) ) for s in 1:length(name) ]
