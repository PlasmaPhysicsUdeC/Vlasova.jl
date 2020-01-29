using Vlasova, Statistics

# BEHAVIOR
@vlasova begin
    data_path = "data/1d_test"
    NUM_THREADS = Threads.nthreads()
    integrator = verlet_velocity
end

# Final conditions
dt = 1e-1
final_time = 500

# Space nodes
box = Box( Nx = 128,
           Nv = 512,
           Lx = 5pi,
           vmin = -6.0,
           vmax = 8.0
           )

# Species declaration
electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = cosine_perturbation1d(box, modes = 1, amplitudes = 5e-2) ⊗ maxwellian1d(box)
                   )


# Make plasma
plasma = Plasma( box = box,
                 species = [ electrons ]
                 )

# Return nothing when this file is included
nothing
