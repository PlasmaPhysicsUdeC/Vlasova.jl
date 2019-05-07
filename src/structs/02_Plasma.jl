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
    species::Array{Specie}
    box::Box
    number_of_species::Int64
    specie_axis::Base.OneTo{Int64}
    
    # Constructors
    # Calculate derived quantities from the array of species
    Plasma(s::Union{T, Array{T}} where T <: Specie,
           p::Box ) = begin
               isarray = typeof(s) <: Array
               isarray ? nothing : (s = [s])
               new(s,
                   p,
                   size(s, 1),
                   Base.OneTo(size(s, 1))
                   )
           end
end
