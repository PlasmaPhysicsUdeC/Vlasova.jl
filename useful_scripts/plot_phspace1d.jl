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

last_it_saved = h5read("electrons.h5", "last_iteration_saved")[1]
times_saved = h5read("electrons.h5", "times_saved")

instant_axis = 1:length(times_saved)

# Portion of phspace to plot
idx = findall( 0 .< box.x[1] .< Lx[1] )
idv = findall( 3 .< box.v[1] .< 4 )

f0 = initial_distribution( box )[idx, idv]

for i in instant_axis
    # Break loop if no posterior times have been saved
    ( times_saved[i] > (last_it_saved - 1) * dt) ? continue : nothing
    
    # Read from file
    fid = h5open("electrons.h5", "r")
    f = fid["distribution"][:,:, i][idx,idv]
    close(fid)

    f = f #.- f0
    
    clf()
    imshow(f', aspect = "auto", origin = "lower",
           extent = [ box.x[1][idx[1]], box.x[1][idx[end]], box.v[1][idv[1]], box.v[1][idv[end]] ],
           interpolation = "bicubic", cmap = cmap )

    colorbar()
    
    my_savefig("at-$(times_saved[i])")
    
    # Done
    println("t = $(times_saved[i]) complete!")
end
