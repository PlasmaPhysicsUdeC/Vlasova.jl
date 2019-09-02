# DEPRECATED

## Simple 1-d simulation

One simple example of a simulation would be

```
using Vlasova

#=    Define parameters	=#

# Name to save the data
simulation_name = "simple_example"

# Space nodes
Nx = 256
Nv = 512

# Space lengths
Lx = 5pi
vMin = -6.0                  # Minimum velocities
vMax = 6.0                   # Maximum velocities

# Final conditions
dt = 5e-2
final_time = 100             # In electron plasma periods

initial_distribution = ones(Nx) * maxwellian1d( box.v[1] ) # maxwellian1d is exported by Vlasova

box = Box(simulation_name, Nx, Nv, Lx, vMin, vMax)
electrons = Specie("electrons", -1.0, 1.0, 1.0, ones(Nx)*maxwellian1d( box ) )

plasma = Plasma( electrons, box )

integrate!(plasma, Nt, dt)
```
