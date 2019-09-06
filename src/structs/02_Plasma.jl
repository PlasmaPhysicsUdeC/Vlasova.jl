"""
```julia
Plasma(;species::Union{Specie, Array{Specie}}, box::Box)
```

Take `species` and `box` and link them into a single container called `Plasma`,
which represents a physical plasma and contains all of the information required
to perform a simulation.

# Notes
* A `Plasma` is a container that points to the `species` and `box`. This means that
  they share the same memory address, and if the plasma is integrated, the distributions in `species`
  will change accordingly.
* To create a `Plasma`, you can provide a single `Specie` or an `Array{Specie}`.
  In the case that a single `Specie` is provided, it will be transformed into an `Array` internally.

# Examples

```jldoctest; setup = :(using Vlasova)
julia> box = Box(Nx = 512,
                 Nv = 1024,
                 Lx = 2pi,
                 vmin = -6,
                 vmax = 6);

julia> electrons = Specie(name = "electrons",
                          charge = -1,
                          mass = 1,
                          temperature = 1,
                          distribution = ones(box.Nx) ⊗ maxwellian1d(box));

julia> plasma = Plasma(species = [electrons],
                       box = box);
```
"""
struct Plasma
    species::Array{Specie, 1}
    box::Box
    number_of_species::Int64
    specie_axis::Base.OneTo{Int64}

    # Constructors
    # Calculate derived quantities from the array of species
    Plasma(args...
           ;species::Union{T, Array{T}} where T <: Specie,
           box::Box ) = begin
               isarray = typeof(species) <: Array
               isarray ? nothing : (species = [species])
               sdims = div.([length( size(species[i].distribution) ) for i in 1:length(species)], 2)
               for i in 2:length(species)
                   @assert sdims[i] == sdims[1] "The sizes of the distribution functions of the species do not match."
               end
               bdim = box.number_of_dims
               sdim = sdims[1]
               bd = "($bdim + $bdim)"
               sd = "($sdim + $sdim)"
               @assert bdim == sdim "The species are $sd-dimensional, but the box provided is $bd-dimensional."
               @assert args == () "Box should not be called using non-keyword arguments."
               new(species,
                   box,
                   size(species, 1),
                   Base.OneTo(size(species, 1))
                   )
           end
end

function Base.display(p::Plasma)
    d = "
    ---
    $(p.box.number_of_dims)-dimensional Vlasova Plasma.
    ---
    "
    print(d)
end
