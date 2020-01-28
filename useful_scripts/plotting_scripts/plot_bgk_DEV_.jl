using Vlasova, PyPlot

box = Box(name = "bgk_creation",
          Nx = 512,
          Nv = 1024,
          Lx = 2pi/0.35,
          vmin = -6,
          vmax = 8);


@time bgk, W = bgk1d(box, amplitude = 0.3, wavenumber = 0.35, vphi = 3.321836);

fignum = 0; close("all")

# z = 0 Profile
figure(fignum += 1)

plot( box.v[1], bgk[1, :] )
ylabel("\$ f(x = 0, v_x) \$")
xlabel("velocity, \$ v_x \$")
grid(true)
axis([2, 5.0, -1e-3, 0.025])


# Phase space filled contour
figure(fignum += 1)

## Filter velocities
idv = findall( 2 .< box.v[1] .< 4.6 );

# Intensity range of interest
vmin = 0.0
vmax = 7e-3

imshow(bgk[:, idv]', aspect = "auto",
                     origin = "lower",
                     extent = [box.x[1][1], box.x[1][end], box.v[1][idv[1]], box.v[1][idv[end]]],
                     vmin = -0.0, vmax = 7e-3, interpolation = "none", cmap = ColorMap("jet"));
colorbar()

xlabel("position, \$ x \$")
ylabel("velocity, \$ v_x \$")

# Distribution as a function of W
figure(fignum += 1)
plot(W[:], bgk[:])
ylabel("\$ f(W) \$")
xlabel("single particle energy, \$ W \$")
grid(true)
