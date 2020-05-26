using Vlasova
using Statistics

# Final conditions
dt = 1e-1
final_time = 6000

# BEHAVIOR

@vlasova begin
    continue_from_backup = true
    data_path = "data/2d_bgk_highres"
    # save_distribution_times = 0:500:Main.final_time |> collect
    integrator = verlet_velocity
    checkpoint_percent = 2
end

# Inject code: Save reduced distribution function
bcode = quote
    @eval using HDF5, Statistics

    # temporaries to calculate momentum density
    tpx = Array{Float64}(undef, plasma.box.Nx..., plasma.box.Nv[1])
    tpy = Array{Float64}(undef, plasma.box.Nx..., plasma.box.Nv[2])

    #fid = h5open(datasaver.path*"/f_reduced.h5", "w")
    ## reduced distribution function
    #fid["fxy"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
    #fid["fxv"] = Array{Float64}(undef, plasma.box.Nx[1], plasma.box.Nv[1], Nt)
    #fid["fyv"] = Array{Float64}(undef, plasma.box.Nx[2], plasma.box.Nv[2], Nt)
    #fid["fvv"] = Array{Float64}(undef, plasma.box.Nv..., Nt)
    ## momentum density
    #fid["px"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
    #fid["py"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
    #close( fid )
end

icode = quote
    tpx .= reducedims(mean, plasma.species[1].distribution, dims = (4))
    Threads.@threads for i in 1:plasma.box.Nv[1]
        @views @. tpx[plasma.box.space_axes..., i] *= plasma.box.v[1][i]
    end

    tpy .= reducedims(mean, plasma.species[1].distribution, dims = (3))
    Threads.@threads for i in 1:plasma.box.Nv[2]
        @views @. tpy[plasma.box.space_axes..., i] *= plasma.box.v[2][i]
    end

    fid = h5open(datasaver.path*"/f_reduced.h5", "r+")
    fid["fxy"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (3, 4))
    fid["fxv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (2, 4))
    fid["fyv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (1, 3))
    fid["fvv"][:, :, t] = reducedims(mean, plasma.species[1].distribution, dims = (1, 2))

    fid["px"][:, :, t] = reducedims(mean, tpx, dims = (3))
    fid["py"][:, :, t] = reducedims(mean, tpy, dims = (3))
    close(fid)
end

Vlasova.inject_to_integrator(before_loop = bcode,
                             inside_loop = icode,
                             after_loop = quote end )


# Enable debugging
# Vlasova.enable_debugging()

# Space nodes
box = Box( Nx = (128, 128),
           Nv = (512, 128),
           Lx = (2pi/0.35, 200pi),
           vmin = (-6.0, -6.0),
           vmax = (8.0, 6.0)
           )

# # Species declaration
# f0 = bgk1d( box, amplitude = 0.3, wavenumber = 0.35, vphi = 3.321836) # ω = 1.16264
# f0 = ones(box.Nx[2]) ⊗ f0
# f0 = permutedims( f0, (2, 1, 3) )

# pert = 1e-15 * rand(box.Nx...)
# pert .-= mean(pert)

# for i in box.space_indices
#     f0[i, :] *= 1 + pert[i]
# end
# f0 = f0 ⊗ maxwellian1d( box.v[2] )


electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = Array{Float64}(undef, box.Nx..., box.Nv...)
                   )


# Make plasma
plasma = Plasma( box = box,
                 species = [ electrons ]
                 )

# Return nothing when this file is included
nothing
