using PyCall, PyPlot
using Vlasova, HDF5
ioff()

include("parameters.jl")

percent = 40
folder = "animation"

# ------------ Filter some data ------------
## Temporally
plot_every = 1

## And spatially
xstep = 3
ystep = 5
# ------------------------------------------

# Make folder to save the plots
run(`mkdir -p $folder`)

np = pyimport("numpy")
# Make space and time axes
x, y = box.x
X, Y = np.meshgrid(x, y)

# Construct time axis
Nt = round(Int, final_time/dt)
last_it = div(percent * Nt, 100 )
time_axis = Array( 1:last_it+1 )*dt

iteration_axis = [ 1 : plot_every : (last_it + 1) ]
times_to_plot = time_axis[ iteration_axis... ]

# Loading calculated quantities
fid = h5open("shared_data.h5", "r")
chargedensity  = fid["chargedensity"][:, :, iteration_axis... ]
close(fid)

Ex, Ey = get_electric_field(chargedensity, box)

# Transpose data
Ex = Ex[1:xstep:end, 1:ystep:end, :]
Ey = Ey[1:xstep:end, 1:ystep:end, :]
chargedensty = chargedensity

println("Generating images:")
let Ex = Ex, Ey = Ey, chargedensity = chargedensity
    xx = x[1:xstep:end]
    yy = y[1:ystep:end]

    extent = [x[1], x[end], y[1], y[end]]
    vmin = -0.005
    vmax =  0.005
    
    for it in 1:length(times_to_plot)
        clf()
        # Arrows
        quiver(xx, yy, Ex[:,:,it]', Ey[:,:,it]', pivot = "mid")

        # Background
        imshow(chargedensity[:,:,it]', aspect = "auto", origin = "lower",
               extent = extent, vmin = vmin, vmax = vmax)#, interpolation = "bilinear")
        
        colorbar()
        
        title("t = $(round(time_axis[it], digits = 1))")
        xlabel("x")
        ylabel("y")
        tight_layout()
        savefig("$folder/$it.png", dpi = 300)

        # Show progress
        print("\r $(round(100it/length(times_to_plot), digits = 3))% done!")
    end
end

println("Generating video:")
run(`ffmpeg -r 15 -framerate 15 -start_number 1 -i $folder/%1d.png $(folder).mp4`)
println("done!")

println("\n\n\nImages will NOT be erased automatically")
