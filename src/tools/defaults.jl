"""
```julia
continue_from_backup = false
```

Specify if Vlasova should look for an unfinished simulation in the current path, and complete it.

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
continue_from_backup = false

"""
```julia
FFTW_flags = FFTW.ESTIMATE
```

Flags used by [`FFTW`](@ref https://github.com/JuliaMath/FFTW.jl) on the Fourier Transform plans inside Vlasova.

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
FFTW_flags = FFTW.ESTIMATE

"""
```julia
NUM_THREADS = 1
```

Number of threads to use in the FFT's inside Vlasova.

# Notes
* This value is passed to [`FFTW`](@ref https://github.com/JuliaMath/FFTW.jl).
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
NUM_THREADS = 1

"""
```julia
integrator = verlet_velocity
```

Default [`symplectic integrator`](@ref https://en.wikipedia.org/wiki/Symplectic_integrator) to be used by Vlasova.

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
integrator = verlet_velocity

"""
```julia
velocity_filtering = true
```

Define if Vlasova should use anisotropic filtering in the transformed-velocity space.

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
velocity_filtering = true

"""
```julia
external_potential(box, time) = 0.0
```

Define some external potential acting on the system. It must be a function of `box::Vlasova.Box` and `time::Float64`

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
external_potential = get_zero

"""
```julia
save_distribution_times = Float64[]
```

Define instants to save the distribution function(s).

# Notes
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
save_distribution_times = Float64[]

"""
```julia
checkpoint_percent = 100
```

Specify some percentage of the simulation to repeatedly back up the data.

# Notes
* When `checkpoint_percent == 100` (default), the simulation will not be backed up.
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
checkpoint_percent = 100


"""
```julia
save_path = "/"
```

Path to save the results of the simulation.

# Notes
* When set to "/" (default value), the results of the simulation will not be saved to the disk.
* This value must be changed from inside Vlasova. If you want to change it, use [`@vlasova`](@ref).
"""
save_path = "/"
