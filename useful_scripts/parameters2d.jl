using Vlasova
using Statistics

# Final conditions
dt = 1e-1
final_time = 500

# BEHAVIOR
@vlasova begin
    data_path = "data/BGK_from_bigger_beam"
    NUM_THREADS = 4
    save_distribution_times = 0:50:500 |> collect
    integrator = verlet_velocity
    checkpoint_percent = 1
end

# Space nodes
box = Box( Nx = (128, 128),
           Nv = (512, 64),
           Lx = (40pi, 200pi),
           vmin = (-6.0, -6.0),
           vmax = (8.0, 6.0)
           )

# Species declaration

let
    fvx = bump_on_tail1d(box.v[1], vtb = 0.5, vdb = 4.5, nc = 0.8, nb = 0.2)
    fvy = maxwellian1d( box.v[2] )
    pertx = 1e-15 * rand(box.Nx[1]) # noise
    perty = 1e-15 * rand(box.Nx[2]) # noise
    pertx .-= mean(pertx)
    perty .-= mean(perty)

    global f0 =  (1 .+ pertx) ⊗ (1 .+ perty) ⊗ fvx ⊗ fvy
end

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
