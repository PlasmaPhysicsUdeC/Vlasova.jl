using Vlasova, Plots

box = Box(Nx = 512,
          Nv = 1024,
          Lx = 2pi/0.35,
          vmin = -6,
          vmax = 8);


@time if !(@isdefined reload_bgk) || reload_bgk
    bgk, W = bgk1d(box, amplitude = 0.3, wavenumber = 0.35, vphi = 3.32)#1836);
    reload_bgk = false
end

# z = 0 Profile
p = plot( box.v[1], bgk[1, :],
          ylabel = "\$ f(x = 0, v_x) \$",
          xlabel = "\$ v_x \$",
          xlims = (2.0, 5.0),
          ylims = (-1e-3, 0.025),
          grid = :true
          )

display(p)

# sleep(10)

## Filter velocities
idv = findall( 2 .< box.v[1] .< 4.6 );

# Intensity range of interest
vmin = 0.0
vmax = 7e-3

p = plot(box.x[1], box.v[1][idv], bgk[:, idv]',
         xlabel = "\$ x \$",
         ylabel = "\$ v_x \$",
         clim = (vmin, vmax),
         colorbar = :true,
         seriestype = :heatmap
         )
display(p)

# Distribution as a function of W
p = plot(W[:], bgk[:],
         ylabel = "\$ f(W) \$",
         xlabel = "\$ W \$",
         grid = :true
         )

display(p)
