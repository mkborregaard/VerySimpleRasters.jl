struct VerySimpleRaster{T}
    mat::Matrix{T}
    nodata::T
    xs::StepRangeLen{Float64}
    ys::StepRangeLen{Float64}
    projection::String
    filename::String
end

VerySimpleRaster(mat::AbstractMatrix{T}, nodata, xll, yll, cellx, celly = cellx, proj = "", fname = "") where T =
    VerySimpleRaster{T}(mat, convert(T, nodata), range(xll, step = cellx, length = size(mat,1)), range(yll, step = celly, length = size(mat,2)), proj, abspath(fname))

bbox(vsr::VerySimpleRaster) = (
    xmin = minimum(vsr.xs),
    xmax = maximum(vsr.xs) + step(vsr.xs),
    ymin = minimum(vsr.ys),
    ymax = maximum(vsr.ys) + step(vsr.ys))

Base.show(io::IO, vsr::VerySimpleRaster{T}) where T = println(io,
"""
VerySimpleRaster{$T} with $(length(vsr.ys)) rows and $(length(vsr.xs)) columns
Resolution: $(round(step(vsr.xs), sigdigits = 5)), $(round(step(vsr.ys), sigdigits = 5))
Extent: $(round(bbox(vsr).xmin, sigdigits = 5)) to $(round(bbox(vsr).xmax, sigdigits = 5)) longitude, $(round(bbox(vsr).ymin, sigdigits = 5)) to $(round(bbox(vsr).ymax, sigdigits = 5)) latitude
Projection: $(vsr.projection)
Source file: $(vsr.filename)""")
