using Vlasova

# Name to save the data
simulation_name = "tst_1d"

# Space nodes
Nx = 256
Nv = 512

# Space lengths
Lx = 5pi
vMin = -6.0       # Minimum velocities
vMax = 8.0        # Maximum velocities

# Final conditions
dt = 1e-1
final_time = 50            # In electron plasma periods

# Save the distribution function ( final_time will always be saved )
save_distribution_times = 0:50:final_time |> collect

# Backup simulation after some % accomplished
checkpoint_percent = 10

# Multi-threding
num_threads = 2

# Multi-specie support
name = ["electrons"]#, "protons"]
charge = [-1.0, 1.0]
mass = [1.0, 1836.15267389]
temperature = [1.0, 1.0]
perturbed = [true, false]

# Choose a Vlasova integrator
integrator = ChinA

# You can comment this function if you are not going to use it
# function external_potential(time, box)

#     cutoff_time = 500
#     cutoff_delay = 20

#     k = 0.35
#     vphi = 3.488
#     potential = @. 0.001 * cos( k * ( box.x[1] - vphi * time ))

#     return @. potential * adiabatic_cutoff( time; cutoff_time = cutoff_time, cutoff_delay = cutoff_delay )
# end

# Starting conditions
function perturbate!(distribution::Array{Float64}, box::Box)
    # Perturbation parameters
    mode = 1
    amplitude = 5e-2

    # Perturbation of the form ( Ax*cos(kx*x) + Ay*cos(ky*y) ... )
    k_mode = @. 2pi * (mode / box.Lx) * box.x
    perturbation = zeros( box.Nx )
    for d in box.dim_axis, i in box.space_indices
        perturbation[i] += amplitude[d] * cos( k_mode[d][ i[d] ] )
    end

    # Ensure that perturbation preserves mass
    perturbation .= perturbation .-  mean(perturbation) .+ 1

    # Apply perturbation
    for j in box.velocity_indices
        for i in box.space_indices
            distribution[i, j] .*= perturbation[i]
        end
    end
    return 0;
end

function initial_distribution(box::Box; perturbate::Bool = false)

    distribution = Array{Float64}(undef, box.N )

    M = Vlasova.maxwellian1d( box.v... )

    for j in CartesianIndices(box.Nv)
        distribution[box.space_axes..., j] .= M[j]
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

# Return nothing when this file is included
nothing
