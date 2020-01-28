# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, Plots, ProgressMeter
import HDF5, FFTW

dirname = "vids"
mkpath( dirname )

period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

# Load files only if this script is open for the first time in this session
@time if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    total_kinetic_energy = HDF5.h5read("shared_data.h5", "total_kinetic_energy")[1:Nt]
    chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]

    # Other quantities
    time_axis = Array([i*dt for i=0:Nt-1]);

    reload_files = false
end

# Wavenumber and frequency
k = rfft_wavevector( box.x )
ω = wavevector( time_axis )

# Reduce quantities
time = time_axis[1:period:Nt]
ρ = chargedensity[:, :, 1:period:Nt]
T = total_kinetic_energy[1:period:Nt]

# Electrostatic energy due to the electric field
V = get_electrostatic_energy(box, ρ)

# Fourier modes of the electrostatic energy density
ρk = FFTW.rfft(ρ, (1,2)) .|> abs2 .|> log

# Plots start from here
timeran = 1:8:Nt
p = Progress(size(timeran, 1), 1) # Progressbar

println("Start loop")

kyind = 1:div(length(k[2]), 2)
anim = @animate for i in timeran
    p1 = plot(time, V,
              xlabel = "\$ t\\omega_{pe} \$",
              ylabel = "\$ V \\frac{ m_e^2 v_{te}^4}{e^2 } \$",
              xlim = (time_axis[1], time_axis[end]),
              ytickfontsize = 6,
              yscale = :log10,
              legend = :false);

    vline!(p1, [time_axis[i]], line = (:red, 1))

    p2 = heatmap( box.x[1], box.x[2], ρ[:,:, i]',
                  title = "\$ \\rho \$",
                  xlabel = "\$ x \\lambda_{D}^{-1} \$",
                  ylabel = "\$ y \\lambda_{D}^{-1} \$",
                  cbar = :bottom,
                  clim = (-0.035, 0.035),
                  c = :plasma)

    p3 = heatmap( k[1], k[2][kyind], ρk[:, kyind, i]',
                  title = "\$ \\log |\\hat \\rho |^2\$",
                  xlabel = "\$ k_x \\lambda_{D} \$",
                  ylabel = "\$ k_y \\lambda_{D} \$",
                  cbar = :bottom,
                  clim = (-15, 5),
                  c = :thermal )

    plot(p1, p2, p3,# p4,
         titlefontsize = 10,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b  c ]
         )

    next!(p)
end

mp4(anim, "$dirname/modes.mp4", fps = 24)
