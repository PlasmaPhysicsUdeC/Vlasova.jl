Any simulation will consist of a set of initial conditions (the **plasma** initial state) and a determined amount of time to integrate it, which will be accomplished through little time steps.

To perform a simulation is as simple as writing

```julia
integrate!(plasma, final_time, dt)
```
where `plasma` holds all the information of the plasma at the initial instant, `final_time` is the amount of time desired to integrate the plasma, normalized to plasma periods of the species used as reference, and `dt` is the time step used to integrate.

## How is a plasma constructed?

In Vlasova, a plasma is an element of type `Plasma`, which itself relies on two structs:

* [`Box`](@ref):
  Represents the phase space occupied by the plasma and the way it is divided to study it.
  Contains all the information regarding the space/velocity extention of the plasma, the number of nodes, etc.

* [`Specie`](@ref):
  Represents one charged species in the plasma. Each species must have a `name`, `charge`, `mass`, `temperature`, and `distribution`.


## Defining the phase space

First, it is necessary to provide the phase space. For example, in the 1-dimensional case one could do

```julia
box = Box(Nx = 512,
          Nv = 2048,
          Lx = 5pi,
          vmin = -6,
          vmax = 6)
```

Where `Nx` and `Nv` are the number of space and velocity nodes, `Lx` is the space extension of the plasma, and `vmin` and `vmax` are the minimum and maximum velocities, normalized to the thermal velocity of each specie, respectively.


## Defining the species that form the plasma

Having the space for the plasma defined, the next step is to define the species that compose the plasma. It is very common to consider only electrons, so let's try that!

```julia
electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = ones(box.Nx) ⊗ maxwellian1d(box) )
```
where [`maxwellian1d`](@ref) is one of the [distributions provided](@ref Distributions) by this same package.

## Making up the plasma

Finally, the plasma is constructed using the information about the space and the species as

```julia
plasma = Plasma(species = [ electrons ],
                box = box)
```

Some important points to keep in mind are:
* The Vlasova integrator will ensure quasineutrality by adding the necessary charge background. In the present case, for example, a neutralizing fixed ion background will be added, and only the electron's dynamics will be considered.
* The plasma struct is just a container *pointing* to the species, meaning that any change performed over `plasma` will be reflected on the species.

## Integrating a plasma

As stated earlier, the plasma may be evolved in time just by typing

```julia
integrate!(plasma, final_time, dt)
```

which will evolve `plasma` until time `final_time` using steps `dt`.

However, there are many orders which can be passed to Vlasova to change the way it works. For example, one may want to save the results to disk, or use checkpoints to backup the data every some percent accomplished of the integration. All this features may be found in the section [Default behavior](@ref).
