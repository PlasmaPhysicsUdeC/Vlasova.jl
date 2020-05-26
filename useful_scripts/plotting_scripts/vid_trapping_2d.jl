# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, Plots, ProgressMeter
using HDF5, FFTW, LinearAlgebra, Strided

function change_reference_frame!(ρ, box, vphi)
    kx, = rfft_wavevector( box.x )

    base_evop = @strided cis.( kx * vphi * dt )        # Evolution operator for one time step
    evop = ones(Complex{Float64}, size(base_evop) )

    FFTW.set_num_threads( Threads.nthreads() )

    plan = plan_rfft(ρ[box.space_axes..., 1], 1)
    @views rhok = plan * ρ[box.space_axes..., 1]

    @inbounds for i in 1:size(ρ)[end]
        @views mul!(rhok, plan, ρ[box.space_axes..., i] )
        @strided rhok .*= evop
        @views ldiv!(ρ[box.space_axes..., i], plan, rhok)
        evop .*= base_evop
    end

    return nothing
end

dirname = "vids"
mkpath( dirname )

period = 1
#period = floor(Int, dt^(-1)); # Uncomment to plot every 1 plasma period

# Load files only if this script is open for the first time in this session
@time if (@isdefined reload_files) ? reload_files : true
    # Load parameters
    include("parameters.jl")

    # Load files
    Nt = HDF5.h5read("electrons.h5", "last_iteration_saved")[1]
    total_kinetic_energy = HDF5.h5read("shared_data.h5", "total_kinetic_energy")[1:Nt]
    chargedensity = HDF5.h5read("shared_data.h5", "chargedensity")[:, :, 1:Nt]

    Lv = prod(box.vmax[i] - box.vmin[i] for i in box.dim_axis)

    fxv= Lv * HDF5.h5read("f_reduced.h5", "fxv")[:, :, 1:Nt]

    # Other quantities
    time_axis = Array([i*dt for i=0:Nt-1]);

    reload_files = false
end

# Reduce quantities
time = time_axis[1:period:Nt]
ρ = chargedensity[:, :, 1:period:Nt]
fxv = fxv[:, :, 1:period:Nt]
T = total_kinetic_energy[1:period:Nt]

# # Change reference frame
# vphi = 3.324 # 3.321836
# change_reference_frame!(ρ, box, vphi)
# change_reference_frame!(px, box, vphi)
# change_reference_frame!(py, box, vphi)

# Electrostatic energy due to the electric field
V = get_electrostatic_energy(box, ρ)

# Plots start from here
timeran = 1:80:Nt # 1:8:Nt
pbar = Progress(size(timeran, 1), 1) # Progressbar

println("Start loop")

anim = @animate for i in timeran
    p1 = plot(time, V,
              xlabel = "\$ t\\omega_{pe} \$",
              ylabel = "\$ V \\frac{ m_e^2 v_{te}^4}{e^2 } \$",
              xlim = (time_axis[1], time_axis[end]),
              ytickfontsize = 6,
              yscale = :log10,
              legend = :false);

    vline!(p1, [time_axis[i]], line = (:red, 1))

    p2 = heatmap( box.x[1], box.x[2], ρ[:,:, i]',
                  title = "\$ \\rho \$",
                  xlabel = "\$ x \\lambda_{D}^{-1} \$",
                  ylabel = "\$ y \\lambda_{D}^{-1} \$",
                  cbar = :bottom,
                  clim = (-0.035, 0.035),
                  c = :plasma)

    p3 = heatmap( box.x[1], box.v[1], fxv[:,:, i]',
                  title = "\$ f_{xv} \$",
                  xlabel = "\$ x \\lambda_{D}^{-1} \$",
                  ylabel = "\$ v v_{te}^{-1} \$",
                  cbar = :bottom,
                  c = :plasma)

    plot(p1, p2, p3,
         titlefontsize = 10,
         dpi = 300,
         layout = @layout [ a{0.2h} ; b c ]
         )

    next!(pbar)
end

mp4(anim, "$dirname/trapping.mp4", fps = 6)
