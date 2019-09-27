using Plots

# start from here
"""
# Notes
* Remember to transpose your input! (or use `transpose = true`).
* The desired output is only produced using the GR backend.
"""
function marginalplot(x, y, a; seriestype = :heatmap, marginalcolor = :black, dpi = 300, kwargs...)

    # reduce along both dimensions
    ax = dropdims(sum(a, dims = 2), dims = 2)
    ay = dropdims(sum(a, dims = 1), dims = 1)

    # upper plot
    u = plot(x, ax,
             xlim = (x[1], x[end]),
             color = marginalcolor,
             legend = :false,
             showaxis = false,
             ticks = nothing)

    # right plot
    r = plot(ay, y,
             ylim = (y[1], y[end]),
             color = marginalcolor,
             legend = :false,
             showaxis = false,
             ticks = nothing)

    # center plot
    c = plot(x, y, a';
             seriestype = seriestype,
             colorbar = false,
             kwargs...)

    return plot(u, c, r, r,
                dpi = dpi,
                layout = @layout [ a{0.1h, 0.9w}  _ ; b{0.9w, 0.9h} c{0.1w} ]
                )
end
