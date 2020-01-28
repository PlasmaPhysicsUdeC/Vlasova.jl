using Vlasova, Statistics

# BEHAVIOR
@vlasova begin
    data_path = "data/1d_BGK_from_beam"
    NUM_THREADS = 2
    integrator = verlet_velocity
end

# Final conditions
dt = 1e-1
final_time = 500

# Space nodes
box = Box( Nx = 128,
           Nv = 512,
           Lx = 40pi,
           vmin = -6.0,
           vmax = 8.0
           )

# Species declaration
let
    fx = 1e-15 * rand( box.Nx[1] )
    fx .-= mean(fx)
    global f0 = ( 1 .+ fx ) ⊗ bump_on_tail1d( box.v[1], vtb = 0.5, vdb = 4.5, nc = 0.9, nb = 0.1 )
end
electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = f0 )


# Make plasma
plasma = Plasma( box = box,
                 species = [ electrons ]
                 )

# Return nothing when this file is included
nothing
