struct VerySimpleRaster{T}
    mat::Matrix{T}
    nodata::T
    xs::StepRangeLen{Float32}
    ys::StepRangeLen{Float32}
    proj::String
    filename::String
end

VerySimpleRaster(mat::AbstractMatrix{T}, nodata, xll, yll, cellx, celly = cellx, proj = "", fname = "") where T =
    VerySimpleRaster{T}(mat, convert(T, nodata), range(xll, step = cellx, length = size(mat,1)), range(yll, step = celly, length = size(mat,2)), proj, abspath(fname))

Base.show(io::IO, vsr::VerySimpleRaster{T}) where T = println(io,
"""
VerySimpleRaster{$T} with $(length(vsr.ys)) rows and $(length(vsr.xs)) columns
Resolution: $(round(vsr.xs.step, sigdigits = 5)), $(round(vsr.ys.step, sigdigits = 5))
Extent: $(round(first(vsr.xs), sigdigits = 5)) to $(round(last(vsr.xs), sigdigits = 5)) longitude, $(round(first(vsr.ys), sigdigits = 5)) to $(round(last(vsr.ys), sigdigits = 5)) latitude
Projection: \"$(vsr.proj)\" """)
Source file: $(vsr.filename)""")
