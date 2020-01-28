# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova
using Plots, ProgressMeter, HDF5

include("../../plotting_scripts/plotting_tools.jl")

dirname = "vids"
mkpath( dirname )

color = :plasma

# Load files only if this script is open for the first time in this session
if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]

    reload_files = false
end

# Time axis
time_axis = collect(0:Nt-1)*dt

# Electrostatic energy
EE = get_electrostatic_energy(box, chargedensity)

# Electrostatic potential
ϕ = get_electrostatic_potential(box, chargedensity)

# Loop
timeran = 30000:50:50000               # Time range to plot

p = Progress(size(timeran, 1), 1) # Progressbar

EEplot = plot(time_axis, EE,
              xlabel = "\$ t\\omega_{pe} \$",
              ylabel = "\$ V \$",
              xlims = (time_axis[1], time_axis[end]),
              ytickfontsize = 6,
              yscale = :log10,
              legend = :false);

anim = @animate for i in timeran

    p1 = vline!(deepcopy(EEplot), [time_axis[i]], line = (:red, 1))

    p2 = contourf(box.x[1], box.x[2], chargedensity[:,:,i]',
                  title = "\$ \\rho \$",
                  xlabel = "\$ x \$",
                  ylabel = "\$ y \$",
                  #colorbar = :true,
                  #clim = (0.004, 0.006),
                  color = color );

    p3 = contourf( box.x[1], box.x[2], ϕ[:, :, i]',
                   title = "\$ \\phi \$",
                   xlabel = "\$ x \$",
                   ylabel = "\$ y \$",
                   #colorbar = :true,
                   #clim = (0.004, 0.006),
                   color = color );

    plot(p1, p2, p3,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b c ]
         )

    next!(p)
end


mp4(anim, "$dirname/density_contour.mp4", fps = 24)
