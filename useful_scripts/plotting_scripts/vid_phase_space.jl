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
    fxy = HDF5.h5read("f_reduced.h5", "fxy")[:,:, 1:Nt]
    fxv = HDF5.h5read("f_reduced.h5", "fxv")[:,:, 1:Nt]
    fyv = HDF5.h5read("f_reduced.h5", "fyv")[:,:, 1:Nt]
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
              ylabel = "\$ V \$",
              xlims = (time_axis[1], time_axis[end]),
              ytickfontsize = 6,
              yscale = :log10,
              legend = :false);

anim = @animate for i in timeran

    p1 = vline!(deepcopy(EEplot), [time_axis[i]], line = (:red, 1))

    p2 = marginalplot(box.x[1], box.x[2], fxy[:,:, i]',
                      xlabel = "\$ x \$",
                      ylabel = "\$ y \$",
                      colorbar = :true,
                      #clim = (0.004, 0.006),
                      color = color );

    p3 = marginalplot(box.x[1], box.v[1], fxv[:,:, i]',
                      xlabel = "\$ x \$",
                      ylabel = "\$ v_x \$",
                      colorbar = :true,
                      #clim = (0.0, 0.005),
                      color = color );

    p4 = marginalplot(box.x[2], box.v[2], fyv[:,:, i]',
                      xlabel = "\$ y \$",
                      ylabel = "\$ v_y \$",
                      colorbar = :true,
                      clim = (0.0, 0.001),
                      color = color );

    p5 = marginalplot(box.v[1], box.v[2], fvv[:,:, i]',
                      xlabel = "\$ v_x \$",
                      ylabel = "\$ v_y \$",
                      colorbar = :true,
                      xlims = (box.v[1][1], box.v[1][end]),
                      ylims = (box.v[2][1], box.v[2][end]),
                      #clim = (0.0, 0.02),
                      color = color );

    plot(p1, p2, p3, p4, p5,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b c ; d e ]
         )

    next!(p)
end


mp4(anim, "$dirname/phase_space.mp4", fps = 24)
