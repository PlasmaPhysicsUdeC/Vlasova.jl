using Vlasova

# Name to save the data
simulation_name = "tst_ABACABA"

# Space nodes
Nx = 256
Nv = 512

# Space lengths
Lx = 5pi
vMin = -6.0                  # Minimum velocities
vMax = 6.0                   # Maximum velocities

# Final conditions
dt = 0.1
final_time = 100             # In electron plasma periods

# Save the distribution function ( final_time will always be saved )
save_distribution_times = [50, 100]

# Multi-threding
num_threads = 2

# Multi-specie support
name = ["electrons"]#, "protons"]
charge = [-1.0, 1.0]
mass = [1.0, 1836.15267389]
temperature = [1.0, 1.0]
perturbed = [true, false]

# Choose a Vlasova integrator

xi = 0.006938106540706989
lambda = 0.2470939580390842
teta = 0.08935804763220157

integrator = VlasovaIntegrator("ABACABA", [teta, lambda, 0.5 - teta, 1 - 2*lambda , 0.5 - teta, lambda, teta],
                               gradient_coefficients = [xi])

# You can comment this function if you are not going to use it
function external_potential(time, box)
    return Vlasova.Zero()
end

# Starting conditions
function perturbate!(distribution::Array{Float64}, box::Box)
    # Perturbation parameters
    mode = 1
    amplitude = 0.05
    
    # Perturbation of the form ( Ax*cos(kx*x) + Ay*cos(ky*y) ... )
    k_mode = @. 2pi * (mode / box.Lx) * box.x
    perturbation = zeros( box.Nx )
    for d in box.dim_axis, i in box.space_indices
        perturbation[i] += amplitude[d] * cos( k_mode[d][ i[d] ] )
    end
    
    # Ensure that perturbation preserves mass
    perturbation = perturbation .-  mean(perturbation) 

    # Apply perturbation
    for i in box.space_indices, j in box.velocity_indices
        distribution[i, j] = ( 1 + perturbation[i] ) * distribution[i, j]
    end
    return 0;
end

function initial_distribution(box::Box; perturbate::Bool = false)

    distribution = ones( box.N )

    for d in box.dim_axis, j in box.velocity_indices, i in box.space_indices
        distribution[i, j] *= Vlasova.maxwellian1d( box.v[d][ j[d] ] )
    end

    # Apply perturbation ?
    perturbate ? perturbate!( distribution, box ) : nothing
    return distribution
end

# Create plasma box
box = Box( name = simulation_name,
           Nx = Nx,
           Nv = Nv,
           Lx = Lx,
           vmin = vMin,
           vmax = vMax );
