"""
    Defines a Vlasovian Plasma, and contains:
        * species: An array of Specie elements
        * box: The parameters of the simulation, of type Box
        * number_of_species: The number of charged species contained in the plasma
        * specie_axis: A tuple of indices, to iterate over each specie.
        
        To be defined, it required to specify species and box:
            julia> parameters = Box( ... );
            julia> plasma = Plasma([Specie(...),
                                    Specie(...)],
                                    parameters );
        *** This struct is dependent on the Specie and Box structs ***
        """
struct Plasma
    species::Array{Specie, 1}
    box::Box
    number_of_species::Int64
    specie_axis::Base.OneTo{Int64}
    
    # Constructors
    # Calculate derived quantities from the array of species
    Plasma(s::Union{T, Array{T}} where T <: Specie,
           b::Box ) = begin
               isarray = typeof(s) <: Array
               isarray ? nothing : (s = [s])
               sdims = div.([length( size(s[i].distribution) ) for i in 1:length(s)], 2)
               for i in 2:length(s)
                   @assert sdims[i] == sdims[1] "The sizes of the distribution functions of the species do not match."
               end
               bdim = b.number_of_dims
               sdim = sdims[1]
               bd = "($bdim + $bdim)"
               sd = "($sdim + $sdim)"
               @assert bdim == sdim "The species are $sd-dimensional, but the box provided is $bd-dimensional."
               new(s,
                   b,
                   size(s, 1),
                   Base.OneTo(size(s, 1))
                   )
           end
end

function Base.display(p::Plasma)
    d = "
    ---
    $(p.box.number_of_dims)-dimensional Vlasova Plasma.
    ---
    
    # TODO: Do display( Specie  ) first
    "
    println(d)
end
