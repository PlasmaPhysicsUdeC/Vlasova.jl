# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, DSP, HDF5, FFTW, Statistics
using Plots

color = :plasma

dirname = "spectrum"
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

# Spectrogram
spect = spectrogram(V, 100, fs = 1 / dt, window = hanning)
freq = 2pi * spect.freq

freqind = Colon()#findall( 0 .<= freq .<= 5 )

heatmap(spect.time, freq[freqind], log10.(spect.power[freqind, :]),
        xlabel = "\$ t \\omega_{pe} \$",
        ylabel = "\$ \\omega \\omega_{pe} \$",
        clim = (-15, 10),
        c = color,
        dpi = 300
        )

savefig("$dirname/spectrogram.png")

# k vs omega
kx, ky = rfft_wavevector( box.x[1] ), rfft_wavevector( box.x[2] )
ω = wavevector( time )

# Fourier modes of the electrostatic energy density
ρx = reducedims(mean, ρ, dims = (2))
ρy = reducedims(mean, ρ, dims = (1))

ρkx = FFTW.rfft(ρx, (1, 2))  .|> abs2 .|> log10
ρky = FFTW.rfft(ρy, (1, 2))  .|> abs2 .|> log10

# ifftshift frequencies and reverse their order because of the sign convention
# 1im ( k*x - ω * t)

ωplot = ifftshift(ω)

#ωind = Colon()
ωind = findall( -5 .< ωplot .< 5 ) # Plot only for some range of freqs

ωplot = ωplot[ωind]
ρkxplot =  reverse(ifftshift(ρkx, 2), dims = 2)[:, ωind]
ρkyplot =  reverse(ifftshift(ρky, 2), dims = 2)[:, ωind]


# rhok(kx, w)
p1 = heatmap(ωplot, kx[2:end-1], ρkxplot[2:end-1, :],
             title = "\$ \\log|\\hat \\rho (k_x, \\omega) |^2 \$",
             xlabel = "\$ \\omega \$",
             ylabel = "\$ k_x \$"
             )

# rhok(kx, w)
p2 = heatmap(ωplot, ky[2:end-1], ρkyplot[2:end-1, :],
             title = "\$ \\log|\\hat \\rho (k_y, \\omega) |^2 \$",
             xlabel = "\$ \\omega \$",
             ylabel = "\$ k_y \$"
             )

plot(p1, p2,
     c = color,
     dpi = 300
     )

savefig("$dirname/disprel.png")
