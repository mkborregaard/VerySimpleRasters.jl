#----------------------------------------
# aggregate

"""
    aggregate(vsr::VerySimpleRaster, factor::Int, fun::Function [, file])

Aggregates the raster by merging `factor` cells in both directions, using the
aggregation function `fun`. Return value is a new raster. Also supports the
`do` block syntax. Use the optional `file` argument to keep the resulting raster
file in your work directory.
"""
aggregate(fun::Function, vsr, factor, file = "") = aggregate(vsr, factor, fun, file)

function aggregate(vsr::VerySimpleRaster, factor::Int, fun::Function, file = "")
   file == "" && (file = tempname())
   file[end-3:end] ∈ (".gri", "grd") && (file = file[1:end-4])
   io = open(file*".gri", "w")
   newm, newn = ceil.(Int,size(vsr.mat)./factor)
   type = typeof(fun([first(vsr.mat)]))
   innerfun(x) = fun(filter(x->x!=vsr.nodata, x))
   for j in 1:newn-1
      for i in 1:newm-1
        write(io, innerfun(vec(@view vsr.mat[(i-1)*factor .+ (1:factor), (j-1)*factor .+ (1:factor)])))
      end
      write(io, innerfun(vec(@view vsr.mat[(newm-1)*factor + 1:size(vsr.mat,1), (j-1)*factor .+ (1:factor)])))
   end
   for i in 1:newm-1
      write(io, innerfun(vec(@view vsr.mat[(i-1)*factor .+ (1:factor), (newn-1)*factor + 1:size(vsr.mat,2)])))
   end
   write(io, innerfun(vec(@view vsr.mat[(newm-1)*factor + 1:size(vsr.mat, 1), (newn-1)*factor + 1:size(vsr.mat,2)])))
   close(io)
   mat = Mmap.mmap(file*".gri", Matrix{type}, (newm, newn))
   vsr_out = VerySimpleRaster(mat, vsr.nodata, minimum(vsr.xs), minimum(vsr.ys), step(vsr.xs)*factor, step(vsr.ys)*factor, vsr.projection, file*".gri")
   write_grd_header(file*".grd", vsr_out)
   vsr_out
end

#--------------------------------------------
# getindex

function getindex(vsr::VerySimpleRaster, x...)
   ret = vsr.mat[x...]
   ret == vsr.nodata && return missing
   ret
end

function getindex(vsr::VerySimpleRaster, x)
   ret = vsr.mat[x]
   ret == vsr.nodata && return missing
   ret
end

#-------------------------------------------
# crop
coordinate_to_index(vsr, x) = coordinate_to_index(vsr, x...)
function coordinate_to_index(vsr, x, y)
   bb = bbox(vsr)
   (x < bb.xmin || x > bb.xmax || y < bb.ymin || y > bb.ymax) && (throw(BoundsError()))
   indx, indy = findfirst(i-> i>x, vsr.xs), findfirst(i-> i>y, vsr.ys)
   indx - 1, length(vsr.ys) - indy + 1
end

index_to_coordinate(vsr, x) = index_to_coordinate(vsr, x...)
function index_to_coordinate(vsr, x, y)
   (x < 1 || x > size(vsr.mat, 1) || y < 1 || y > size(vsr.mat, 2)) && throw(BoundsError())
   vsr.xs[x]+step(vsr.xs), vsr.ys[length(vsr.ys)-y+1]
end

"""
    crop(vsr::VerySimpleRaster{T}, xmin, xmax, ymin, ymax [, file])
    crop(vsr::VerySimpleRaster{T}, extent [, file])

Crops the raster to a rectangular window. Returns a new raster, use the
`file` argument to keep this in your working directory.
"""
crop(vsr, inds, file = "") = crop(vsr, inds..., file)
function crop(vsr::VerySimpleRaster{T}, xmin, xmax, ymin, ymax, file = "") where T
   x1,y1 = coordinate_to_index(vsr, xmin, ymax)
   x2,y2 = coordinate_to_index(vsr, xmax, ymin)
   newsize = (x2-x1, y2-y1)

   file == "" && (file = tempname())
   file[end-3:end] ∈ (".gri", "grd") && (file = file[1:end-4])
   open(file*".gri", "w") do IO
      write(IO, vsr.mat[x1:x2-1, y1+1:y2])
   end
   mat = Mmap.mmap(file*".gri", Matrix{T}, newsize)
   vsr_out = VerySimpleRaster(mat, vsr.nodata, vsr.xs[x1], vsr.ys[length(vsr.ys)-y2], step(vsr.xs), step(vsr.ys), vsr.projection, file*".gri")
   write_grd_header(file*".grd", vsr_out)
   vsr_out
end

#------------------------------------------
# extract

"""
   extract(vsr::VerySimpleRaster, x, y)
   extract(vsr::VerySimpleRaster, points)

Extracts the value of the raster at the given points.
"""
extract(vsr::VerySimpleRaster, x, y) = extract(vsr, (x, y))
function extract(vsr::VerySimpleRaster, tup)
   x, y = tup[1], tup[2]
   indx, indy = try
      coordinate_to_index(vsr, x, y)
   catch
      return missing
   end
   vsr[indx, indy]
end

extract(vsr::VerySimpleRaster, tup::Tuple{T, T}) where T <: AbstractVector = extract.(Ref(vsr), zip(tup...))
extract(vsr::VerySimpleRaster, ar::AbstractVector) = extract.(Ref(vsr), ar)

#-------------------------------------------
# mask

midpoints(a::AbstractVector) = [mean(a[i-1:i]) for i in eachindex(a)[2:end]]
const Pt{T<:Real} = Tuple{T, T}

function isinside(r::Pt, poly::AbstractVector)
    # An implementation of Hormann-Agathos (2001) Point in Polygon algorithm
    # See: http://www.sciencedirect.com/science/article/pii/S0925772101000128
    # Code segment adapted from PolygonClipping.jl
    c = false
    detq(q1,q2,r) = (q1[1]-r[1])*(q2[2]-r[2])-(q2[1]-r[1])*(q1[2]-r[2])
    for i in eachindex(poly)[2:end]
        q2 = poly[i]
        q1 = poly[i-1]
        if q1 == r
            @warn("point on polygon vertex - returning false")
            return false
        end
        if q2[2] == r[2]
            if q2[1] == r[1]
                @warn("point on polygon vertex - returning false")
                return false
            elseif (q1[2] == r[2]) && ((q2[1] > x) == (q1[1] < r[1]))
                @warn("point on edge - returning false")
                return false
            end
        end
        if (q1[2] < r[2]) != (q2[2] < r[2]) # crossing
            if q1[1] >= r[1]
                if q2[1] > r[1]
                    c = !c
                elseif ((detq(q1,q2,r) > 0) == (q2[2] > q1[2])) # right crossing
                    c = !c
                end
            elseif q2[1] > r[1]
                if ((detq(q1,q2,r) > 0) == (q2[2] > q1[2])) # right crossing
                    c = !c
                end
            end
        end
    end
    return c
end

"""
    mask(raster, polygon [, filename])`

Masks the raster by a polygon. Creates a new raster where points falling outside
the polygon have been replaced by `missing`
"""
function mask(vsr::VerySimpleRaster{T}, poly::AbstractVector, file = "") where T
   file == "" && (file = tempname())
   file[end-3:end] ∈ (".gri", "grd") && (file = file[1:end-4])
   open(file*".gri", "w") do IO
      for j in axes(vsr.mat, 2)
         for i in axes(vsr.mat, 1)
            if isinside(index_to_coordinate(vsr, i, j), poly)
               write(IO, vsr[i, j])
            else
               write(IO, vsr.nodata)
            end
         end
      end
   end
   mat = Mmap.mmap(file*".gri", Matrix{T}, size(vsr.mat))
   vsr_out = VerySimpleRaster(mat, vsr.nodata, first(vsr.xs), first(vsr.ys), step(vsr.xs), step(vsr.ys), vsr.projection, file*".gri")
   write_grd_header(file*".grd", vsr_out)
   vsr_out
end
