using Vlasova

# BEHAVIOR
@vlasova begin
    data_path = "data/valentini_2005"
    NUM_THREADS = 2
    integrator = mclachlan_velocity
    checkpoint_percent = 10
end

# Final conditions
dt = 1e-1
final_time = 200

# Space nodes
box = Box( Nx = 256,
           Nv = 512,
           Lx = 5pi,
           vmin = -6.0,
           vmax = 6.0
           )

# Species declaration
electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = cosine_perturbation1d(box, modes = 1, amplitudes = 5e-2) ⊗ maxwellian1d( box )
                   )


# Make plasma
plasma = Plasma( box = box,
                 species = [ electrons ]
                 )

# Return nothing when this file is included
nothing
