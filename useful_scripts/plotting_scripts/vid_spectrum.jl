# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, HDF5, FFTW, Statistics, ProgressMeter
using Plots

color = :plasma

dirname = "vids"
mkpath( dirname )

# Load files only if this script is open for the first time in this session
@time if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    ρ = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]

    # Time axis
    time = collect(0:Nt-1)*dt

    reload_files = false
end

# Electrostatic energy
V = get_electrostatic_energy(box, ρ)

# k vs omega
kx = rfft_wavevector( box.x[1] )
ky = rfft_wavevector( box.x[2] )

# Reduced charge density
ρx = reducedims(mean, ρ, dims = 2)
ρy = reducedims(mean, ρ, dims = 1)

function get_spectrogram( ρr; window = 50 )
    Nx, Nt = size(ρr)
    wρ = Array{Float64}(undef, div(Nx, 2) + 1, 2*window+1, Nt)

    for i in (window+1):(Nt-window)
        wρ[:, :, i] = rfft(ρr[:, (i-window):(i+window)], (1, 2) ) .|> abs2
    end

    return wρ
end

window = 250

sρx = get_spectrogram( ρx, window = window )
sρy = get_spectrogram( ρy, window = window )


# Prepare quantities to plot
ω = wavevector( Array(1:(2window + 1))*dt )

ωplot = ifftshift( ω )
ωind = findall(-5 .< ωplot .< 5)
ωplot = ωplot[ ωind ]

sρx_plot = ifftshift( sρx, 2 )
sρy_plot = ifftshift( sρy, 2 )

# Invert frequencies to fulfill sign convention: 1im ( k * x - ω * t)
sρx_plot = reverse( sρx_plot, dims = 2)[:, ωind, :] .|> log10
sρy_plot = reverse( sρy_plot, dims = 2)[:, ωind, :] .|> log10

# Electrostatic energy plot
p0 = plot(time, V,
          xlabel = "\$ t\\omega_{pe} \$",
          ylabel = "\$ V \$",
          xlims = (time[1], time[end]),
          ytickfontsize = 6,
          yscale = :log10,
          legend = :false);

timeran = (1+window):5:(Nt-window)

p = Progress(size(timeran, 1))

anim = @animate for i in timeran
    p1 = vline!( deepcopy(p0), [time[i]], line = (:red, 1) )

    vspan!(p1, [time[i] - dt*window, time[i] + dt*window],
           color = :gray,
           alpha = 0.3)

    # Ex
    p2 = heatmap(ωplot, kx[2:end-1], sρx_plot[2:end-1, :, i],
                 title = "\$ \\log| \\hat \\rho (k_x, \\omega) |^2 \$",
                 xlabel = "\$ \\omega \$",
                 ylabel = "\$ k_x \$"
                 )

    p3 = heatmap(ωplot, ky[2:end-1], sρy_plot[2:end-1, :, i],
                 title = "\$ \\log| \\hat \\rho (k_y, \\omega) |^2 \$",
                 xlabel = "\$ \\omega \$",
                 ylabel = "\$ k_y \$"
                 )

    plot(p1, p2, p3,
         c = color,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b c ]
         )

    next!(p)
end

mp4(anim, "$dirname/spectrum.mp4", fps = 24)
