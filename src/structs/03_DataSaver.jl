mutable struct DataSaver
    chargedensity::Array{Float64}
    kinetic_energy::Array{Float64}
    checkpoint_axis::Array{Float64, 1}
    last_iteration_saved::Int64
    last_distribution_saved::Int64
    checkpoints_reached::Int64
    save_distribution_times::Array{Float64, 1}
    path::String
    
    DataSaver(plasma::Plasma,
              Nt::Integer,
              checkpoint_percent::Integer = 10,
              continue_from_backup::Bool = false,
              save_distribution_times::Array{Float64, 1} = Float64[]
              ) = begin

                  # checkpoints
                  @assert 1 <= checkpoint_percent <= 100  "checkpoint_percent must be valued between 1 and 100"
                  @assert isinteger(100/checkpoint_percent) "checkpoint_percent must be a whole divisor of 100"
                  checkpoint_axis = collect(Int, range(1, stop = Nt, length = Int(100/checkpoint_percent) + 1 ))
                  checkpoint_step = checkpoint_axis[2] - checkpoint_axis[1]

                  chargedensity = Array{Float64}(undef, (plasma.box.Nx..., checkpoint_step) )
                  kinetic_energy = Array{Float64}(undef, (checkpoint_step, plasma.number_of_species))
                  last_iteration_saved = 1
                  last_distribution_saved = 0
                  checkpoints_reached = 1
                  
                  path = "data/"*plasma.box.simulation_name*"/"
                  Ndf = length( save_distribution_times )

                  if continue_from_backup
                      # Restore data
                      for s in 1:plasma.number_of_species
                          fid = HDF5.h5open(path*plasma.species[s].name*".h5", "r")
                          if s == 1
                              last_iteration_saved = read( fid["last_iteration_saved"])[1]
                              last_distribution_saved = read( fid["last_distribution_saved"] )[1]
                          end
                          plasma.species[s].distribution .= read( fid["distribution"],
                                                                  (UnitRange.(1, plasma.box.N)..., last_distribution_saved) )
                          HDF5.close(fid)
                      end
                      checkpoints_reached = findfirst( last_iteration_saved .== checkpoint_axis )
                  else
                      # Initialize files
                      # Make folder to save data
                      mkpath(path)
                      
                      # Common file
                      fid = HDF5.h5open(path*"shared_data.h5", "w")
                      fid["chargedensity"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
                      fid["total_kinetic_energy"] = Array{Float64}(undef, Nt)
                      fid["specie_names"] = [plasma.species[s].name for s in plasma.specie_axis]
                      ## Save first instant
                      fid["chargedensity"][UnitRange.(1, plasma.box.Nx)..., 1] = get_density( plasma )
                      fid["total_kinetic_energy"][1] = sum( get_kinetic_energies(plasma) )
                      HDF5.close(fid)

                      # Specie file(s)
                      save_df0 = ( 0.0 in save_distribution_times )
                      for s in 1:plasma.number_of_species
                          fid = HDF5.h5open(path*plasma.species[s].name*".h5", "w")
                          fid["distribution"] = Array{Float64}(undef, plasma.box.N..., Ndf)

                          if save_df0
                              fid["distribution"][UnitRange.(1, plasma.box.N)..., 1] = plasma.species[s].distribution
                          end
                          
                          fid["times_saved"] = save_distribution_times
                          if s == 1
                              fid["last_iteration_saved"] = [1]
                              fid["last_distribution_saved"] = ( save_df0 ? [1] : [0] )
                          end
                          HDF5.close(fid)
                      end
                      save_df0 ? ( last_distribution_saved += 1 ) : nothing
                  end
                  
                  new( chargedensity,
                       kinetic_energy,
                       checkpoint_axis,
                       last_iteration_saved,
                       last_distribution_saved,
                       checkpoints_reached,
                       save_distribution_times,
                       path)
              end
end
