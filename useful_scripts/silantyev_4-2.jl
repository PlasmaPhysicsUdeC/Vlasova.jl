using Vlasova
using FFTW
using Plots; pyplot()


box = Box(Nx = (64, 64),
          Nv = (256, 32),
          Lx = (2pi/0.35, 400pi),
          vmin = (-6, -6),
          vmax = (6, 6)
          );

f01d = bgk1d(box, wavenumber = 0.35, amplitude = 0.2, vphi = 3.3582, dim = 1)[1]
f0 = permutedims( ones(box.Nx[1]) ⊗ f01d, (2, 1, 3)) ⊗ maxwellian1d(box.v[2])

electrons = Specie(name = "electrons",
                   charge = -1,
                   mass = 1,
                   temperature = 1,
                   distribution = f0
                   );

dens = get_density(box, electrons)

Ex, Ey = get_electric_field(box, dens)
Ekx = rfft( Ex, (1, 2) ) .|> abs

p =  Ekx[2, :] / length(Ekx) .+ 1e-16

p1 = scatter(p, yscale = :log10, ylim = [1e-5, 1])
hline!(p1, [1e-3], lw = 2)

EE = get_electrostatic_energy(box, dens)
K = get_kinetic_energy(box, electrons )
TE = EE + K
EE = EE / TE
K = K / TE


p2 = scatter([EE], scale = :log10, label = "EE")
scatter!(p2, [K], label = "K")

plot(p1, p2)

# return
#nothing
