# Workaround for https://github.com/JuliaPlots/Plots.jl/issues/1905
ENV["GKSwstype"]="100"

using Vlasova, HDF5, Plots


dirname = "plots"
mkpath( dirname )

include("parameters.jl")

times_saved = h5read("electrons.h5", "times_saved")
instant_axis = 1:12#:length(times_saved)

fxy = Array{Array{Float64, 2}, 1}(undef, length(instant_axis))
fvv = Array{Array{Float64, 2}, 1}(undef, length(instant_axis))
fxv = Array{Array{Float64, 2}, 1}(undef, length(instant_axis))
fyv = Array{Array{Float64, 2}, 1}(undef, length(instant_axis))
for i in instant_axis
    # Read from file
    fid = h5open("electrons.h5", "r")
    f = fid["distribution"][:,:,:,:, i][:,:,:,:]
    close(fid)

    #f = f .- f0

    # fxy
    fxy[i] = reducedims(sum, f, dims = box.velocity_dims ) / prod( box.Nv )

    # fvv
    fvv[i] = reducedims(sum, f, dims = box.space_dims ) / prod( box.Nx )

    # fxv
    fxv[i] = reducedims(sum, f, dims = (2, 4) ) / ( box.Nx[2] * box.Nv[2] )

    # fyv
    fyv[i] = reducedims(sum, f, dims = (1, 3) ) / ( box.Nx[1] * box.Nv[1] )

    println("t = ", times_saved[i], "/", times_saved[end], " calculated!")
end
f = nothing # Release memory
println("\nNow, plotting...")

# Plots
color = :plasma

begin
    # VV
    pxy = [ heatmap(box.x[1],box.x[2], fxy[i]',
                    title = "\$ \\omega_{pe} t  = $(times_saved[i]) \$",
                    colorbar = :false,
                    xticks = :false,
                    yticks = :false,
                    color = color)
            for i in instant_axis]

    pvv = [ heatmap(box.v[1],box.v[2], fvv[i]',
                    title = "\$ \\omega_{pe} t = $(times_saved[i]) \$",
                    colorbar = :false,
                    xticks = :false,
                    yticks = :false,
                    color = color)
            for i in instant_axis]

    pxv = [ heatmap(box.x[1],box.v[1], fxv[i]',
                    title = "\$ \\omega_{pe} t = $(times_saved[i]) \$",
                    colorbar = :false,
                    xticks = :false,
                    yticks = :false,
                    color = color)
            for i in instant_axis]

    pyv = [ heatmap(box.x[2],box.v[2], fyv[i]',
                    title = "\$ \\omega_{pe} t = $(times_saved[i]) \$",
                    colorbar = :false,
                    xticks = :false,
                    yticks = :false,
                    color = color)
            for i in instant_axis]

    # Put scale and ticks

    plot!(pxy[9],
          xlabel = "\$ x \\lambda_{D}^{-1} \$",
          ylabel = "\$ y \\lambda_{D}^{-1} \$",
          plot_title = "\$ \\rho \$",
          xticks = :true,
          yticks = :true)

    plot!(pvv[9],
          xlabel = "\$ v_x v_{te}^{-1} \$",
          ylabel = "\$ v_y v_{te}^{-1} \$",
          xticks = -6:2:8,
          yticks = -6:2:6)

    plot!(pxv[9],
          xlabel = "\$ x \\lambda_{D}^{-1} \$",
          ylabel = "\$ v_x v_{te}^{-1} \$",
          xticks = :true,
          yticks = :true)

    plot!(pyv[9],
          xlabel = "\$ y \\lambda_{D}^{-1} \$",
          ylabel = "\$ v_y v_{te}^{-1} \$",
          xticks = 0:200:600,
          yticks = :true)

    # Plot and save figs

    plot(pxy...,
         dpi = 500,
         titlefont = 8)

    savefig("$dirname/xy.png")

    plot(pvv...,
         clim = (0.0, 0.002),
         cbar = :true,
         dpi = 500,
         titlefont = 8)

    savefig("$dirname/vv.png")

    plot(pxv...,
         clim = (0.0, 0.001),
         dpi = 500,
         titlefont = 8)

    savefig("$dirname/xv.png")

    plot(pyv...,
         #clim = (0.0, 0.001),
                  cbar = :true,
         dpi = 500,
         titlefont = 8
         )

    savefig("$dirname/yv.png")
end
