using Plots

# start from here
"""
# Notes
* Remember to transpose your input! (or use `transpose = true`).
* The function has only been tested using the GR backend.
"""
function marginalplot(x, y, a; seriestype = :heatmap, marginalcolor = :black, dpi = 300, kwargs...)

    # reduce along both dimensions ( a is transposed )
    ax = dropdims(sum(a, dims = 1), dims = 1)
    ay = dropdims(sum(a, dims = 2), dims = 2)

    uxlim = (x[1], x[end])
    uylim = let min = minimum(ax), max = maximum(ax)
        isapprox(min, max) ? (min - 1e-3, max + 1e-3) : (min, max)
    end
    rxlim = let min = minimum(ay), max = maximum(ay)
        isapprox(min, max) ? (min - 1e-3, max + 1e-3) : (min, max)
    end
    rylim = (y[1], y[end])

    # upper plot
    u = plot(x, ax,
             xlim = uxlim,
             ylim = (minimum(ax), maximum(ax)),
             color = marginalcolor,
             legend = :false,
             showaxis = false,
             ticks = nothing)

    # right plot
    r = plot(ay, y,
             xlim = (minimum(ay), maximum(ay)),
             ylim = rylim,
             color = marginalcolor,
             legend = :false,
             showaxis = false,
             ticks = nothing)

    # center plot
    c = plot(x, y, a;
             seriestype = seriestype,
             colorbar = false,
             kwargs...)

    return plot(u, c, r, r,
                dpi = dpi,
                margins = -5Plots.PlotMeasures.mm ,
                layout = @layout [ a{0.1h, 0.9w}  _ ; b{0.9w, 0.9h} c{0.1w} ]
                )
end
