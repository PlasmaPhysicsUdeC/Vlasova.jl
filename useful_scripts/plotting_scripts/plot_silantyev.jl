# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, Plots, Statistics
import HDF5, FFTW

dirname = "plots"
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

# Reduce quantities
time = time_axis[1:period:Nt]
ρ = chargedensity[:, :, 1:period:Nt]
T = total_kinetic_energy[1:period:Nt]

#electrostatic_energy = get_electrostatic_energy( chargedensity, box );
Efield = get_electric_field(box, ρ)

# Wavenumbers and frequency
k = rfft_wavevector( box.x )
ω = wavevector( time )

# ====================== Plots start from here
mode_x = 1
modes_y = Array( 0:7 )

# Julia arrays start at 1
mode_x += 1
@. modes_y = modes_y + 1

Ekx = FFTW.rfft(Efield[1], (1, 2))[mode_x, : , :] .|> abs2

p = plot()
for mode in modes_y
    Ekxmode = Ekx[ mode, :]
    #gamma = find_exponential_growth(time, Ekxmode, interval = [1400, 1900])

    plot!(p, time, Ekxmode, label = "\$ k_y = $(round(k[2][mode], digits = 3)) \$", lw = 2)#, \\gamma = $(round(gamma, digits = 6))")
end
plot!(p,
      dpi = 300,
      legend = :bottomright,
      grid = :true,
      yscale = :log10,
      xlabel = "\$ t\\omega_{pe} \$",
      ylabel = "\$ | E_x (k_x = $(round(k[1][mode_x], digits = 3)), k_y, t)|^2 \$",
      xlim = (time_axis[1], time_axis[end])
      )
savefig(p, "$dirname/harmonics_growth.png")


# Energy density
edens = reducedims(sum, abs2.(Efield[1]), dims = 1 )*box.dx[1] / box.Lx[1]

p = heatmap(time, box.x[2], log10.(edens),
            title = "\$ \\log \\left( L_x^{-1}\\int |E_x|^2 dx \\right) \$",
            titlefontsize = 12,
            xlabel = "\$ t\\omega_{pe} \$",
            ylabel = "\$ y \\lambda_{De} \$",
            cbar = :true,
            dpi = 300,
            clim = (-15, maximum(log10.(edens)))
            )
savefig("$dirname/energy_density.png")
