using Vlasova, HDF5, PyPlot
ioff()

# Define new function to save figures. It sets tight margins, saves the plot and closes the current figure.
function my_savefig(name::String)
    tight_layout()
    savefig(folder*"/"*name*".png", dpi = 300)
    close( plt.gcf().number )
    clf()
    return nothing;
end

cmap = ColorMap("jet")      # Color map for the filled contour plots

include("parameters.jl")


folder = "plot_phspace"
run(`mkdir -p $folder`)

times_saved = h5read("electrons.h5", "times_saved")

instant_axis = 1:length(times_saved)


run(`mkdir -p $folder/xy`)
run(`mkdir -p $folder/xv`)
run(`mkdir -p $folder/yv`)
run(`mkdir -p $folder/vv`)

f0 = initial_distribution( box )

for i in instant_axis
    # Read from file
    fid = h5open("electrons.h5", "r")
    f = fid["distribution"][:,:,:,:, i][:,:,:,:]
    close(fid)

    f = f .- f0

    # Integrate f along space dims
    fxy = prod(box.dv) * reducedims(sum, f, dims = box.velocity_dims )
    fvv = prod(box.dx) * reducedims(sum, f, dims = box.space_dims )

    fxv = box.dx[2] * box.dv[2] * reducedims(sum, f, dims = (2, 4) )
    fyv = box.dx[1] * box.dv[1] * reducedims(sum, f, dims = (1, 3) )
    
    # Release memory
    f = nothing 
    
    # fxy
    clf()
    imshow(fxy', aspect = "auto", origin = "lower",
           extent = [ box.x[1][1], box.x[1][end], box.x[2][1], box.v[2][end] ],
           interpolation = "none", cmap = cmap )

    colorbar()
    
    title("t = $(times_saved[i])")
    xlabel("Position, x")
    ylabel("Position, y")
    
    my_savefig("xy/xy-at-$(times_saved[i])")

    # fvv
    clf()
    imshow(fvv', aspect = "auto", origin = "lower",
           extent = [ box.v[1][1], box.v[1][end], box.v[2][1], box.v[2][end] ],
           interpolation = "none", cmap = cmap )

    colorbar()

    title("t = $(times_saved[i])")
    xlabel("Velocity, v_x")
    ylabel("Velocity, v_y")
    
    my_savefig("vv/vv-at-$(times_saved[i])")

    # fxv
    clf()
    imshow(fxv', aspect = "auto", origin = "lower",
           extent = [ box.x[1][1], box.x[1][end], box.v[1][1], box.v[1][end] ],
           interpolation = "none", cmap = cmap )

    colorbar()
    
    title("t = $(times_saved[i])")
    xlabel("Position, x")
    ylabel("Velocity, v_x")
    
    my_savefig("xv/xv-at-$(times_saved[i])")

    # fyv
    clf()
    imshow(fyv', aspect = "auto", origin = "lower",
           extent = [ box.x[2][1], box.x[2][end], box.v[2][1], box.v[2][end] ],
           interpolation = "none", cmap = cmap )
    
    colorbar()

    title("t = $(times_saved[i])")
    xlabel("Position, y")
    ylabel("Velocity, v_y")
    
    my_savefig("yv/yv-at-$(times_saved[i])")
   
    # Done
    println("t = $(times_saved[i]) complete!")
end
