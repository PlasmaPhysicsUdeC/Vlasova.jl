"""
    Advects a plasma in space efficiently using the information stored in an
    object of type SpaceAdvection.
"""
function (space_advection::SpaceAdvection)(plasma::Plasma; advection_number::Integer = 1)
    for s in 1:plasma.number_of_species
        LinearAlgebra.mul!(space_advection.transformed_DF,
                           space_advection.plan, plasma.species[s].distribution)
        # Here happens the magic!
        _space_advection!( space_advection.transformed_DF,
                           space_advection.shift[advection_number][s],
                           space_advection.N2p1...)
        
        LinearAlgebra.ldiv!(plasma.species[s].distribution,
                            space_advection.plan, space_advection.transformed_DF)
    end
    return 0;
end

# 1D - Fortran
function _space_advection!(transformed_DF, shift, Nx::Int32, Nvx::Int32)
    ccall((:space_advection1d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Array{Complex{Float64}, 2}}, Ptr{Array{Complex{Float64}, 2}}),
          Ref(Nx), Ref(Nvx), shift, transformed_DF)
end

# 2D - Fortran
function _space_advection!(transformed_DF, shift, Nx::Int32, Ny::Int32, Nvx::Int32, Nvy::Int32)
    ccall((:space_advection2d_, "advections.so"), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Array{Complex{Float64}, 4}}, Ptr{Array{Complex{Float64}, 4}}),
          Ref(Nx), Ref(Ny), Ref(Nvx), Ref(Nvy), shift, transformed_DF)
    return 0;
end
