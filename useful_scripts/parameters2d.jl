using Vlasova
using Statistics

# Final conditions
dt = 1e-1
final_time = 500

# BEHAVIOR
@vlasova begin
    data_path = "data/2d_BGK_from_pert"
    NUM_THREADS = 4
    integrator = verlet_velocity
    checkpoint_percent = 1
end

# Inject code: Save reduced distribution function
bcode = quote
    @eval using HDF5, Statistics
    fid = h5open(datasaver.path*"/f_reduced.h5", "w")

    fid["fxy"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
    fid["fxv"] = Array{Float64}(undef, plasma.box.Nx[1], plasma.box.Nv[1], Nt)
    fid["fyv"] = Array{Float64}(undef, plasma.box.Nx[2], plasma.box.Nv[2], Nt)
    fid["fvv"] = Array{Float64}(undef, plasma.box.Nv..., Nt)

    close( fid )
end

icode = quote
    fid = h5open(datasaver.path*"/f_reduced.h5", "r+")
    fid["fxy"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (3, 4))
    fid["fxv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (2, 4))
    fid["fyv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (1, 3))
    fid["fvv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (1, 2))
    close(fid)
end

Vlasova.inject_to_integrator(before_loop = bcode,
                             inside_loop = icode,
                             after_loop = quote end)

# Space nodes
box = Box( Nx = (128, 128),
           Nv = (512, 64),
           Lx = (40pi, 200pi),
           vmin = (-6.0, -6.0),
           vmax = (8.0, 6.0)
           )

# Species declaration

let
    fvx = maxwellian1d( box.v[1] )
    fvy = maxwellian1d( box.v[2] )
    pertx = cosine_perturbation1d(box, modes = 1, amplitudes = 5e-2)
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
