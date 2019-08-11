using Vlasova

# Name to save the data
simulation_name = "silantyev_non-adiabatic"

# Space nodes
Nx = (64, 64)
Nv = (256, 32)

# Space lengths
Lx = (2pi/0.35, 800pi)
vMin = (-6, -6)       # Minimum velocities
vMax = ( 8,  6)       # Maximum velocities

# Final conditions
dt = 1e-1
final_time = 6000            # In electron plasma periods

# Save the distribution function ( final_time will always be saved )
save_distribution_times = 0:100:6000 |> collect

# Backup simulation after some % accomplished
checkpoint_percent = 5

# Multi-threding
num_threads = 8

# Multi-specie support
name = ["electrons"]#, "protons"]
charge = [-1.0, 1.0]
mass = [1.0, 1836.15267389]
temperature = [1.0, 1.0]
perturbed = [true, false]

# Choose a Vlasova integrator
integrator = verlet_velocity

function Vlasova.anisotropic_filter(box::Box) # Silantyev hyper-viscosity

    D16 = 1e-25
    
    wavevector = rfft_wavevector( box.v )
    Nv2p1, fourier_indices = get_rfft_dims( box.v )
    
    filter = ones( Nv2p1 )
    for i in fourier_indices
        filter[i] = exp.( - dt * D16 * wavevector[1][i[1]]^16 )  # Only damp modes for vx
    end

    return filter
end

# You can comment this function if you are not going to use it
function external_potential(time, box) # Integrated while the src code is not fixed #TODO
    cutoff_time = 110
    cutoff_delay = 0

    k = 0.35
    vphi = 3.488
    w = k * vphi
    
    pot = @. dt * 0.01 * (-1 / w) * cos( k * box.x[1] - w  * (time + dt/4) ) * ones( box.Nx[2] )' # Only right for velocity Verlet
    
    return @. pot * adiabatic_cutoff( time; cutoff_time = cutoff_time, cutoff_delay = cutoff_delay )
end

# Starting conditions
function perturbate!(distribution::Array{Float64}, box::Box)
    # Perturbation parameters
    mode = (0, 0)
    amplitude = (0, 0)
    
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
    
    M2d = Vlasova.maxwellian2d( box.v... )

    for j in box.velocity_indices
        for i in box.space_indices
            distribution[i, j] .= M2d[j]
        end
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
