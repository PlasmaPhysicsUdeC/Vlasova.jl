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

# To rasterize a plot: gca()[:set_rasterization_zorder](-1)
#
# Example:
# gca()[:set_rasterization_zorder](-1)
# contourf(x, v, f'[:,:], zorder = -2)

include("parameters.jl")

period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

# Loading calculated quantities
Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, 1:Nt]
total_kinetic_energy = HDF5.h5read("shared_data.h5", "total_kinetic_energy")[1:Nt]

time_axis = Array([i*dt for i=0:Nt-1]);

k = rfft_wavevector( box.x[1] )
w = wavevector( time_axis )

electrostatic_energy = get_electrostatic_energy( chargedensity, box );

# Reduce quantities
time_axis = time_axis[1:period:Nt]
chargedensity = chargedensity[:, 1:period:Nt]
electrostatic_energy = electrostatic_energy[1:period:Nt]
total_kinetic_energy = total_kinetic_energy[1:period:Nt]


######################## Plots ########################
fig_number = 0

# Electric field energy
fig_number += 1; figure(fig_number)
let
    semilogy(time_axis, electrostatic_energy, linestyle="-", marker="");
    ylabel("Electrostatic energy");
    xlabel("Time, \$\\omega_{pe}^{-1}\$");
    grid(true)
    xlim(time_axis[1], time_axis[end])
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
    xlim(time_axis[1], time_axis[end])
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
    xlim(time_axis[1], time_axis[end])
    my_savefig("totalEnergy")
end

# # Power spectrum of the electrostatic energy
fig_number += 1; figure(fig_number)

let
    
    energy_density_kw = get_power_per_mode(chargedensity, box)
    
    modes_to_plot = 1:4
    for mode in modes_to_plot
        semilogy(time_axis, energy_density_kw[mode+1,:], linestyle="-", marker="", linewidth=0.5, label=string("mode ",mode) );
    end
    
    leg = legend(loc="lower right");
    ylabel("Electrostatic energy per spatial mode");
    xlabel("Time, \$\\omega_{pe}^{-1}\$");

    # change linewidth in legend
    for line in leg.get_lines()
        line.set_linewidth(2)
    end
    grid(true)
    xlim(time_axis[1], time_axis[end])
    my_savefig("power_per_mode")
end

# # Dispersion relation

fig_number += 1; figure(fig_number);
disprel = get_dispersion_relation(chargedensity, box)

# let
#     spaceWaveVector = rfftWaveVector( space ); kMax = maximum(spaceWaveVector);
#     frequency = rfftWaveVector( time_axis );

#     transformedElectricField =  FFTW.rfft(electricField, [2, 1]); # rfft in time and fft in space

#     toPlot = log.( abs.( transformedElectricField ) );
#     toPlot = FFTW.fftshift(toPlot, 1);

#     imshow( toPlot'[:,:] , aspect = "auto", origin="lower", extent = [-kMax, kMax, frequency[1], frequency[end]], interpolation = "none", zorder = -2, cmap = cmap)
#     colorbar();
#     xlabel("Wavenumber, \$ k\\lambda_{De}^{-1}\$")
#     ylabel("Frequency, \$ \\omega \\omega_{pe} \$")
#     my_savefig("dispersionRelation")
# end

# #close("all")

println("All graphics have been plotted.")
