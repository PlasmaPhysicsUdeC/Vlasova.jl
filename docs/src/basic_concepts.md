## Basic concepts

Any simulation will consist of a set of initial conditions (the *plasma* initial state) and a determined ammount of time to integrate it, which will be accomplished through little time steps.

To perform a simulation is as simple as writing

```
integrate(plasma, final_time, dt)
```
where *plasma* holds all the information of the plasma at the initial instant, *final_time* is the ammount of time desired to integrate the plasma, normalized to electron plasma periods, and *dt* is the *time step* used to integrate.

But, how is the *plasma* entity constructed?

In Vlasova, a plasma relies on two structs:

* Box:
  Represents the physical space occupied by the plasma and the way it is divided to study it.
  Contains all the information regarding the space/velocity extention of the plasma, the number of nodes, etc.

* Specie:
  Represents a specie of charged particle in the plasma. Each specie must have a name, charge, mass, temperature, and distribution, and a plasma may be formed by many species.


## Creating a plasma

First, it is necessary to provide the physical space:

```
box = Box(simulation_name, Nx, Nv, Lx, vmin, vmax)
```

Where *simulation_name* is a String, used to identify the simulation later, *Nx* and  *Nv* are the number of space and velocity nodes, *Lx* is the space extention of the plasma, normalized in electrod Debye lengths, and *vmin* and *vmax* are the minimum and maximum velocities, normalized to the thermal velocity of each specie.

Please note that in the 1-d case, *Nx*, *Nv*, *Lx*, *vmin* and *vmax* are usual numbers, but for the multi-dimensional case, they are Tuples or Arrays.


Then, a specie or array of species must be created. An example may be

```
species = [ Specie("some_specie",  some_charge,  some_mass,  some_temperature,  some_initial_distribution ),
	    Specie("other_specie", other_charge, other_mass, other_temperature, other_initial_distribution) ]
```

where all of the specie properties must be normalized to the electron respective magnitudes.

And finally, the plasma may be built as

```
plasma = Plasma(species, box)
```

Test:
This is ``\LaTeX`` code: `` f(x) = \frac{3}{2} x^2 ``
```math
\LaTeX
```
