@recipe function f(val::VerySimpleRaster)
    seriestype --> :heatmap
    aspect_ratio --> 1
    grid --> false
    replace(val.mat, missing => NaN)
end
