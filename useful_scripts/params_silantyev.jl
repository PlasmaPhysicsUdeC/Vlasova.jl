using Vlasova
using Statistics

# Final conditions
dt = 1e-1
final_time = 10

# BEHAVIOR
@vlasova begin
    data_path = "data/2d_silantyev"
    save_distribution_times = 0:500:Main.final_time |> collect
    integrator = verlet_velocity
    checkpoint_percent = 2
end

# Space nodes
box = Box( Nx = (64, 64),
           Nv = (256, 32),
           Lx = (2pi/0.35, 200pi),
           vmin = (-8.0, -6.0),
           vmax = (6.0, 6.0)
           )

# Species declaration

f0 = bgk1d( box, amplitude = 0.3, wavenumber = 0.35, vphi = 3.324,
            wave_frame = true)
f0 = ones(box.Nx[2]) ⊗ f0
f0 = permutedims( f0, (2, 1, 3) )

pert = 1e-15 * rand(box.Nx...)
pert .-= mean(pert)

for i in box.space_indices
    f0[i, :] *= 1 + pert[i]
end
f0 = f0 ⊗ maxwellian1d( box.v[2] )


electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = f0
                   )


# Make plasma
plasma = Plasma( box = box,
                 species = [ electrons ]
                 )

# Return nothing when this file is included
nothing
