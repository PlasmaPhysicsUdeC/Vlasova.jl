mutable struct DataSaver
    chargedensity::Array{Float64}
    kinetic_energy::Array{Float64}
    checkpoint_axis::Array{Float64, 1}
    save_distribution_axis::Array{Float64, 1}
    last_iteration_saved::Int64
    last_distribution_saved::Int64
    checkpoints_reached::Int64
    path::String
    save_data::Bool

    DataSaver(plasma::Plasma,
              Nt::Integer, dt::Float64,
              path::String,
              checkpoint_percent::Integer,
              continue_from_backup::Bool,
              save_distribution_times::Array{N, 1} where N <: Union{Float64, Int64}
              ) = begin

                  # Save data
                  save_data = !isempty(data_path) ||
                      !isempty(save_distribution_times) ||
                      ( checkpoint_percent < 100 )

                  if isempty( data_path )
                      @assert checkpoint_percent == 100 "The variable data_path must be specified if you want to use checkpoints."
                      @assert !continue_from_backup "The variable data_path must be specified if you want to continue from a backup."
                      @assert isempty( save_distribution_times ) "The variable data_path must be specified if you want to save distributions."
                  end

                  # Save distribution axis
                  save_distribution_times = sort( Float64.( save_distribution_times ) )
                  # Always save the last distribution function ( slot needed for backing up )
                  append!(save_distribution_times, (Nt-1)*dt)
                  ## All saving times must be multiple of dt
                  @assert all(isinteger, save_distribution_times ./ dt) "Not all times to save the distribution function are multiples of dt."
                  @assert all(save_distribution_times .<= (Nt-1)*dt ) "The times to save the distribution function can not be larger than final_time."

                  save_distribution_axis = round.(Int, save_distribution_times ./ dt) .+ 1
                  Ndf = size( save_distribution_axis, 1 )

                  # Checkpoint axis
                  @assert 1 <= checkpoint_percent <= 100  "checkpoint_percent must be valued between 1 and 100."
                  @assert isinteger(100/checkpoint_percent) "checkpoint_percent must be a whole divisor of 100."
                  checkpoint_axis = collect(Int, range(1, stop = Nt, length = Int(100/checkpoint_percent) + 1 ))
                  checkpoint_step = checkpoint_axis[2] - checkpoint_axis[1]

                  chargedensity = Array{Float64}(undef, (plasma.box.Nx..., checkpoint_step) )
                  kinetic_energy = Array{Float64}(undef, (checkpoint_step, plasma.number_of_species))
                  last_iteration_saved = 1
                  checkpoints_reached = 1
                  last_distribution_saved = 0

                  if continue_from_backup
                      @assert isfile(joinpath(path, plasma.species[1].name*".h5")) "There is no valid backup in the specified folder."
                      # Restore data
                      for s in plasma.specie_axis
                          fid = HDF5.h5open( joinpath(path, plasma.species[s].name*".h5"), "r")
                          if s == 1
                              last_iteration_saved = fid["last_iteration_saved"][1][1]
                              last_distribution_saved = fid["last_distribution_saved"][1][1]
                          end
                          plasma.species[s].distribution .= dropdims( fid["distribution"][UnitRange.(1, plasma.box.N)..., Ndf],
                                                                      dims = 2 * plasma.box.number_of_dims + 1 )
                          HDF5.close(fid)
                          @assert last_iteration_saved < Nt "The simulation saved in the folder '$path' has already ended successfully."
                      end
                      checkpoints_reached = findfirst( last_iteration_saved .== checkpoint_axis )

                  elseif save_data
                      # Initialize files
                      mkpath(path)

                      # Common file
                      fid = HDF5.h5open( joinpath(path, "shared_data.h5"), "w")
                      fid["chargedensity"] = Array{Float64}(undef, plasma.box.Nx..., Nt)
                      fid["total_kinetic_energy"] = Array{Float64}(undef, Nt)
                      fid["specie_names"] = [plasma.species[s].name for s in plasma.specie_axis]
                      ## Save first instant
                      fid["chargedensity"][UnitRange.(1, plasma.box.Nx)..., 1] = get_density( plasma )
                      fid["total_kinetic_energy"][1] = get_kinetic_energy(plasma)
                      HDF5.close(fid)

                      # Specie file(s)
                      save_df0 = ( 0.0 in save_distribution_times )
                      save_df0 ? ( last_distribution_saved += 1 ) : nothing

                      for s in plasma.specie_axis
                          fid = HDF5.h5open(joinpath(path, plasma.species[s].name*".h5"), "w")
                          fid["distribution"] = Array{Float64}(undef, plasma.box.N..., Ndf)

                          if save_df0
                              fid["distribution"][plasma.box.distribution_axes..., 1] = plasma.species[s].distribution
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
                       path,
                       save_data)
              end
end
