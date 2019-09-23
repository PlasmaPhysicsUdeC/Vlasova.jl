using Vlasova

# BEHAVIOR
@vlasova begin
    data_path = "data/valentini_2005_injected"
    NUM_THREADS = 2
    integrator = mclachlan_velocity
    checkpoint_percent = 10
end

# Injected Code
bcode = quote                   # Before loop
    @eval using HDF5, Statistics
    fid = h5open(datasaver.path*"/f_reduced.h5", "w")
    fid["fx"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
    fid["fv"] = Array{Float64}(undef, plasma.box.Nv..., Nt)
    fid["fx"][:, 1] = reducedims(mean, plasma.species[1].distribution, dims = 2)
    fid["fv"][:, 1] = reducedims(mean, plasma.species[1].distribution, dims = 1)
    close(fid)
end

icode = quote                   # Inside loop
    fid = h5open(datasaver.path*"/f_reduced.h5", "r+")
    fid["fx"][:, t] = reducedims(mean, plasma.species[1].distribution, dims = 2)
    fid["fv"][:, t] = reducedims(mean, plasma.species[1].distribution, dims = 1)
    close(fid)
end

acode = quote                   # After loop
    println("All done!")
end

Vlasova.inject_to_integrator(before_loop = bcode,
                             inside_loop = icode,
                             after_loop = acode )

# Debugging
Vlasova.enable_debugging()

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
