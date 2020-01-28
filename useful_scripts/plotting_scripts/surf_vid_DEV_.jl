# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova
using Plots, ProgressMeter, HDF5

color = :plasma

# Load files only if this script is open for the first time in this session
if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    fvv = HDF5.h5read("f_reduced.h5", "fvv")[:,:, 1:Nt]

    chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]

    reload_files = false
end

# Time axis
time_axis = collect(0:Nt-1)*dt

# Electrostatic energy
EE = get_electrostatic_energy(box, chargedensity)

# Loop
timeran = 1:1:Nt               # Time range to plot

p = Progress(size(timeran, 1), 1) # Progressbar

EEplot = plot(time_axis, EE,
              xlabel = "\$ t\\omega_{pe} \$",
              ylabel = "\$ E^2 \$",
              xlims = (time_axis[1], time_axis[end]),
              ytickfontsize = 6,
              yscale = :log10,
              legend = :false);


anim = @animate for i in timeran

    p1 = vline!(deepcopy(EEplot), [time_axis[i]], line = (:red, 1))

    p2 = surface(box.v[1], box.v[2], fvv[:,:, i]',
                 xlabel = "v_x",
                 ylabel = "v_y",
                 colorbar = :false,
                 zlim = (0.0, 0.15),
                 #camera = (10, 60),
                 alpha = 0.0,
                 grid = true,
                 color = color );

    p3 = surface(box.v[1], box.v[2], fvv[:,:, i]',
                 xlabel = "v_x",
                 ylabel = "v_y",
                 colorbar = :false,
                 clim = (0.0, 0.03),
                 camera = (40, 50),
                 grid = true,
                 color = color );

    plot(p1, p2, p3,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b c ]
         )

    next!(p)
end


mp4(anim, "velocity_surf.mp4", fps = 16)
