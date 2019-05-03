using Vlasova, PyPlot
import HDF5, FFTW

# ioff()
println("\n ---------------------------------------------------------")
println(" Interactive windows are disabled by default [ ioff() ].")
println("   To enable interactivity use the function: ion()")
println("\n ---------------------------------------------------------")

run(`mkdir -p plots`)           # Folder to save the plots

# Define new function to save figures. It sets tight margins, saves the plot and closes the current figure.
my_savefig(name::String) = begin
    tight_layout()
    savefig(name)
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

# Loading calculated quantities
fid = HDF5.h5open("shared_data.h5", "r")
chargedensity, total_kinetic_energy, specie_names = HDF5.read(fid, "chargedensity", "total_kinetic_energy", "specie_names");
close(fid)

@show hasnan(chargedensity)

#
Nt = size(total_kinetic_energy, 1)
timeAxis = Array([i*dt for i=0:Nt-1]);


period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

electrostatic_energy = get_electrostatic_energy( chargedensity, box );


######################## Plots ########################
fig_number = 0

# # Power spectrum of the electric field
# fig_number += 1; figure(fig_number)

# let # To avoid polluting the global scope
#     spacePowerSpectrum = abs2.( FFTW.rfft( electricField, [1, 2] ) );

#     modesToPlot = 1:4
#     for mode in modesToPlot
#         semilogy(timeAxis[1:period:end], spacePowerSpectrum[mode+1, 1:period:end], linestyle="-", marker="", linewidth=0.5, label=string("mode ",mode));
#     end

#     leg = legend(loc="lower right");
#     ylabel("Electrostatic energy per spatial mode");
#     xlabel("Time, \$\\omega_{pe}^{-1}\$");

#     # change linewidth in legend
#     for line in leg.get_lines()
#         line.set_linewidth(2)
#     end
#     grid(true)
#     my_savefig("plots/electrostaticPowerPerMode.pdf")
# end


# Electric field energy
fig_number += 1; figure(fig_number)
let
    semilogy(timeAxis[1:period:end], electrostatic_energy[1:period:end], linestyle="-", marker="");
    ylabel("Electrostatic energy");
    xlabel("Time, \$\\omega_{pe}^{-1}\$");
    grid(true)
    my_savefig("plots/electrostatic_energy.pdf")
end

# Kinetic energy
fig_number += 1; figure(fig_number)
let
    plot(timeAxis[1:period:end], total_kinetic_energy[1:period:end], linestyle="-", marker="", color="red");
    title("Kinetic energy")
    xlabel("Time, \$\\omega_{pe}^{-1}\$")
    ylabel("Kinetic energy")
    grid(true)
    my_savefig("plots/kineticEnergy.pdf")
end

# Total energy
fig_number += 1; figure(fig_number)
let
    totalEnergy = electrostatic_energy + total_kinetic_energy;

    plot(timeAxis[1:period:end], totalEnergy[1:period:end], linestyle="-", marker="", color="red");
    xlabel("Time, \$\\omega_{pe}^{-1}\$")
    ylabel("Total energy, \$ E_p + E_k \$")
    grid(true)
    my_savefig("plots/totalEnergy.pdf")
end

# # Dispersion relation

# fig_number += 1; figure(fig_number);
# let
#     spaceWaveVector = rfftWaveVector( space ); kMax = maximum(spaceWaveVector);
#     frequency = rfftWaveVector( timeAxis );

#     transformedElectricField =  FFTW.rfft(electricField, [2, 1]); # rfft in time and fft in space

#     toPlot = log.( abs.( transformedElectricField ) );
#     toPlot = FFTW.fftshift(toPlot, 1);

#     imshow( toPlot'[:,:] , aspect = "auto", origin="lower", extent = [-kMax, kMax, frequency[1], frequency[end]], interpolation = "none", zorder = -2, cmap = cmap)
#     colorbar();
#     xlabel("Wavenumber, \$ k\\lambda_{De}^{-1}\$")
#     ylabel("Frequency, \$ \\omega \\omega_{pe} \$")
#     my_savefig("plots/dispersionRelation.pdf")
# end

# #close("all")

# println("All graphics have been plotted.")
