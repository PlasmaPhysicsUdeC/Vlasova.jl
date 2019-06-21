mutable struct DataSaver
    chargedensity::Array{Float64}
    kinetic_energy::Array{Float64}
    checkpoint_axis::Array{Float64, 1}
    save_distribution_axis::Array{Float64, 1}
    last_iteration_saved::Int64
    last_distribution_saved::Int64
    checkpoints_reached::Int64
    path::String
    
    DataSaver(plasma::Plasma,
              tm::TimeManager,
              checkpoint_percent::Integer = 10,
              continue_from_backup::Bool = false,
              save_distribution_times::Array{N, 1} where N <: Union{Float64, Int64} = Float64[]
              ) = begin

                  Nt = floor(Int, tm.final_time / tm.dt + 1)

                  # Save distribution axis
                  save_distribution_times = sort( Float64.( save_distribution_times ) )
                  ## All saving times must be multiple of dt
                  @assert all( isinteger.( round.(save_distribution_times ./ tm.dt, digits = 10) )
                               ) "Not all times to save the distribution function are multiples of dt"
                  @assert all(save_distribution_times .- tm.final_time .<= 0) "The times to save the distribution function can not be larger than final_time"
                  ## Always save the last instant
                  (tm.final_time in save_distribution_times) ? nothing : append!(save_distribution_times, tm.final_time)
                  
                  save_distribution_axis = floor.(Int, save_distribution_times ./ tm.dt .+ 1)
                  Ndf = size( save_distribution_axis, 1 )
                  
                  # Checkpoint axis
                  @assert 1 <= checkpoint_percent <= 100  "checkpoint_percent must be valued between 1 and 100"
                  @assert isinteger(100/checkpoint_percent) "checkpoint_percent must be a whole divisor of 100"
                  checkpoint_axis = collect(Int, range(1, stop = Nt, length = Int(100/checkpoint_percent) + 1 ))
                  checkpoint_step = checkpoint_axis[2] - checkpoint_axis[1]

                  chargedensity = Array{Float64}(undef, (plasma.box.Nx..., checkpoint_step) )
                  kinetic_energy = Array{Float64}(undef, (checkpoint_step, plasma.number_of_species))
                  last_iteration_saved = 1
                  checkpoints_reached = 1
                  last_distribution_saved = 0

                  # Path
                  path = "data/"*plasma.box.simulation_name*"/"

                  if continue_from_backup
                      # Restore data
                      for s in 1:plasma.number_of_species
                          fid = HDF5.h5open(path*plasma.species[s].name*".h5", "r")
                          if s == 1
                              last_iteration_saved = read( fid["last_iteration_saved"])[1]
                              last_distribution_saved = read( fid["last_distribution_saved"] )[1]
                          end
                          plasma.species[s].distribution .= read( fid["distribution"],
                                                                  (UnitRange.(1, plasma.box.N)..., Ndf) )
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
                      save_df0 ? ( last_distribution_saved += 1 ) : nothing
                      
                      for s in 1:plasma.number_of_species
                          fid = HDF5.h5open(path*plasma.species[s].name*".h5", "w")
                          fid["distribution"] = Array{Float64}(undef, plasma.box.N..., Ndf)

                          if save_df0
                              fid["distribution"][UnitRange.(1, plasma.box.N)..., 1] = plasma.species[s].distribution
                          end
                          
                          fid["times_saved"] = save_distribution_times
                          if s == 1
                              fid["last_iteration_saved"] = [ last_iteration_saved ]
                              fid["last_distribution_saved"] = [ last_distribution_saved ]
                          end
                          HDF5.close(fid)
                      end
                  end
                  
                  new( chargedensity,
                       kinetic_energy,
                       checkpoint_axis,
                       save_distribution_axis,
                       last_iteration_saved,
                       last_distribution_saved,
                       checkpoints_reached,
                       path)
              end
end
