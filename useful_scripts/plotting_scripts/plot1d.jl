# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, Plots
import HDF5

dirname = "plots"
mkpath( dirname )

period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

# Load files only if this script is open for the first time in this session
if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    total_kinetic_energy = HDF5.h5read("shared_data.h5", "total_kinetic_energy")[1:Nt]
    chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, 1:Nt]

    # Other quantities
    time_axis = Array([i*dt for i=0:Nt-1]);

    reload_files = false
end

# Reduce quantities
time = time_axis[1:period:Nt]
ρ = chargedensity[:, 1:period:Nt]
T = total_kinetic_energy[1:period:Nt]

#electrostatic_energy = get_electrostatic_energy( chargedensity, box );
Efield = get_electric_field(box, ρ)

# EStatic energy due to the electric field
V = 0.5 * prod(box.dx) * reducedims(sum, abs2.( Efield[1] ), dims = (1))

######################## Plots ########################

# Electric field energy
p = plot(time, V, color = :red, legend = :false)

plot!(p,
      dpi = 300,
      grid = :true,
      yscale = :log10,
      xlim = (time[1], time[end]),
      xlabel = "\$ t\\omega_{pe} \$",
      ylabel = "Electrostatic energy")

savefig(p, "$dirname/electrostatic_energy.png")


# Kinetic energy
p = plot(time, T, color = :red, legend = :false)
plot!(p,
      dpi = 300,
      grid = :true,
      xlim = (time[1], time[end]),
      xlabel = "\$ t\\omega_{pe} \$",
      ylabel = "Kinetic energy")

savefig(p, "$dirname/kinetic_energy.png")

# Total energy
Etot = T + V
Etot = ( Etot .- Etot[1] ) ./ Etot[1]

p = plot(time, Etot, color = :red, legend = :false)
plot!(p,
      dpi = 300,
      grid = :true,
      xlim = (time[1], time[end]),
      xlabel = "\$ t\\omega_{pe} \$",
      ylabel = "\$ \\frac{E(t) - E(0)}{E(0)} \$")

savefig(p, "$dirname/total_energy.png")

println("All graphics have been plotted.")
