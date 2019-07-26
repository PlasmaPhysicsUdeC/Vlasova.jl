using Vlasova, PyPlot
import HDF5, FFTW

ioff()                          # Disable interactive graphics
run(`mkdir -p plots`)           # Folder to save the plots

# Define new function to save figures. It sets tight margins, saves the plot and closes the current figure.
function my_savefig(name::String)
    tight_layout()
    savefig("plots/"*name*".png", dpi = 300)
    close( plt.gcf().number )
    return nothing;
end

cmap = ColorMap("jet")      # Color map for the filled contour plots

@time include("parameters.jl")

period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

# Loading calculated quantities
Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]
total_kinetic_energy = HDF5.h5read("shared_data.h5", "total_kinetic_energy")[1:Nt]

@show hasnan(chargedensity)

time_axis = Array([i*dt for i=0:Nt-1]);

k = rfft_wavevector( box.x[1] )
w = wavevector( time_axis )

# Reduce quantities
time_axis = time_axis[1:period:Nt]
chargedensity = chargedensity[:, :, 1:period:Nt]
total_kinetic_energy = total_kinetic_energy[1:period:Nt]

#electrostatic_energy = get_electrostatic_energy( chargedensity, box );
Ex, Ey = get_electric_field(chargedensity, box)

# EStatic energy due to x and y electric field component
EEx = 0.5*prod(box.dx) * reducedims(sum, abs2.(Ex), dims = (1,2))
EEy = 0.5*prod(box.dx) * reducedims(sum, abs2.(Ey), dims = (1,2))

## Total EStatic energy
electrostatic_energy = EEx .+ EEy

######################## Plots ########################
fig_number = 0

# Electric field energy
fig_number += 1; figure(fig_number)
let
    semilogy(time_axis, electrostatic_energy, "b-", label = "Total");
    semilogy(time_axis, EEx, "r--", label = "Ex");
    semilogy(time_axis, EEy, "g--", label = "Ey");
    
    ylabel("Electrostatic energy");
    xlabel("Time, \$\\omega_{pe}^{-1}\$");
    grid(true)
    legend(loc="best")
    my_savefig("electrostatic_energy")
end

# Kinetic energy
fig_number += 1; figure(fig_number)
let
    plot(time_axis, total_kinetic_energy, linestyle="-", marker="", color="red");
    title("Kinetic energy")
    xlabel("Time, \$\\omega_{pe}^{-1}\$")
    ylabel("Kinetic energy")
    grid(true)
    my_savefig("kineticEnergy")
end

# Total energy
fig_number += 1; figure(fig_number)
let
    totalEnergy = electrostatic_energy + total_kinetic_energy;

    plot(time_axis, totalEnergy, linestyle="-", marker="", color="red");
    xlabel("Time, \$\\omega_{pe}^{-1}\$")
    ylabel("Total energy, \$ E_p + E_k \$")
    grid(true)
    my_savefig("totalEnergy")
end

println("All graphics have been plotted.")
