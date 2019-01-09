struct VerySimpleRaster
    mat::Matrix{Union{Missing, Float64}}
    xs::StepRangeLen{Float64}
    ys::StepRangeLen{Float64}
end

VerySimpleRaster(mat, xll, yll, cell) =
    VerySimpleRaster(mat, StepRangeLen(xll, cell, size(mat,2)), StepRangeLen(yll, cell, size(mat,1)))

Base.show(io::IO, vsr::VerySimpleRaster) = println(io,
"""
VerySimpleRaster{Float64} with $(length(vsr.ys)) rows and $(length(vsr.xs)) columns
Extent: $(round(first(vsr.xs), sigdigits = 6)) to $(round(last(vsr.xs), sigdigits = 6)) longitude, $(round(first(vsr.ys), sigdigits = 6)) to $(round(last(vsr.ys), sigdigits = 6)) latitude""")
