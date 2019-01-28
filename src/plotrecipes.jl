@recipe function f(val::VerySimpleRaster)
    seriestype --> :heatmap
    aspect_ratio --> 1
    grid --> false
#    framestyle --> :none
    val.xs, val.ys, @view replace(val.mat, val.nodata => NaN)'[end:-1:1, :]
end
