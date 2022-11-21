# VerySimpleRasters.jl
Very simple raster format for julia

## NOTE: All functionality of this package has been superceded by the Rasters.jl package - I recommend using that package instead

The VerySimpleRasters package allows you to open on-disk rasters as a memory-mapped file, and do simple data curation and abstraction methods on it. It currently uses the native `.grd` format of R's `raster` package and can import ESRI Ascii grids to that format. Further integration with GDAL is planned, which would allow for compatibility with all raster formats.

The VerySimpleRaster type is a simple wrapper around an mmapped `Array` with some metadata. All operations happen on disk, creating temporary on-disk copies. Providing the optional `filename` argument creates the copy permanently at the given path.

Currently available functions are:
#### For loading and writing files
- `VerySimpleRaster(filename)` loads an R .grd file
- `importASCII(filename)` imports an ESRI Ascii grid to a .grd file and opens it
- `writeraster(filename, raster)` writes the raster as a .grd file

#### Raster operations
- `crop(raster, window [, filename])` crops the raster to a window
- `extract(raster, points)` extracts the value of the raster at points
- `mask(raster, polygon [, filename])` masks the raster by a polygon
- `aggregate(raster, factor, fun [, filename])` aggregates the raster by merging `factor` cells in both directions, using aggregation function `fun`

## Sample analysis
#### Load relevant packages and load a raster of global mean temperature
```julia
using VerySimpleRasters, Plots
default(seriescolor = :topo) #not the most correct choice of color
bio1 = VerySimpleRaster("temperature/bio1.grd")
# VerySimpleRaster{Float32} with 900 rows and 2160 columns
# Resolution: 0.16667, 0.16667
# Extent: -180.0 to 180.0 longitude, -60.0 to 90.0 latitude
# Projection: +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0
# Source file: /Users/michael/bio1.gri

plot(bio1)
```
<img width="886" alt="skaermbillede 2019-01-29 kl 08 25 50" src="https://user-images.githubusercontent.com/8429802/51891456-859de300-239f-11e9-8f1f-f310eb83de91.png">

#### Crop the raster to a window
```julia
samerica = crop(bio1, (-80, -35, -60, 15))
plot(samerica)
```
<img width="474" alt="skaermbillede 2019-01-29 kl 08 25 08" src="https://user-images.githubusercontent.com/8429802/51891468-8c2c5a80-239f-11e9-952c-1014d345e1bc.png">

#### Extract part of the raster by a Shapefile polygon
This is more complicated, as Shapefile does not handle dbf files and easy subsetting. So we currently need to jump through some hoops. This will be improved in the future.
```julia
using Shapefile, DBFTables
shp = open("countries/CNTRY92.shp") do IO
    read(IO, Shapefile.Handle)
end
dbf = open("countries/CNTRY92.dbf") do IO
    DBFTables.read_dbf(IO)
end
brazil = shp.shapes[dbf.NAME .== "Brazil"][1]

plot(samerica)
plot!(brazil, fa = 0, lc = :red)
```
<img width="465" alt="Samerica with Brazil" src="https://user-images.githubusercontent.com/8429802/51890627-cfd19500-239c-11e9-8e2b-2cd78f0e46d5.png">

Current issues: VerySimpleRasters doesn't yet know about Shapefiles. It just knows about polygons, and that they are a vector of Tuples, so first I have to extract the points as a Tuple. Also, this is a multipart polygon, so I define a function to only extract the first.
```julia


function splitpolygon(poly::AbstractVector{T}) where T<:Tuple
    ret = Vector{SubArray{T,1}}()
    current = 1
    while current < length(poly)
        next = findnext(x->x == poly[current], poly, current + 1)
        isnothing(next) && (next = length(poly))
        push!(ret, view(poly, current:next))
        current = next
    end
    ret
end

bp = [(pt.x, pt.y) for pt in brazil.points]
brapolys = splitpolygon(bp)
```
Then I can `mask` it
```julia
brarast = mask(samerica, brapolys[1])
# VerySimpleRaster{Float32} with 450 rows and 270 columns
# Resolution: 0.16667, 0.16667
# Extent: -80.0 to -35.0 longitude, -60.0 to 15.0 latitude
# Projection: +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0
# Source file: /var/folders/b8/rvm2dk495f74088ssh9v2l8r0000gn/T/julianZEtiV.gri
```
#### Extract points from that new raster
First sample some random locations from the bounding box
```julia
minx, maxx = extrema(p[1] for p in bp)
miny, maxy = extrema(p[2] for p in bp)
x,y = minx.+(maxx-minx).*rand(50), miny.+(maxy-miny).*rand(50)
```
Then plot it and extract the value under the points
```julia
plot(brarast)
scatter!(x,y, shape = :+, msc = :red, legend = false)
```
<img width="478" alt="Brazil with points" src="https://user-images.githubusercontent.com/8429802/51890641-d7913980-239c-11e9-88b2-b7eb78d365db.png">

```julia
value = extract(brarast, x, y)
value[1:5]
# 5-element Array{Union{Missing, Float32},1}:
#  269.0f0   
#  259.0f0   
#     missing
#  188.0f0   
#     missing
```

#### Aggregate values to a coarser resolution
Here we use the optional filename argument to keep the resulting raster in the working directory

```julia
using Statistics
bra_coarse = aggregate(brarast, 9, mean, "aggregate_result")
plot(bra_coarse)
```

<img width="473" alt="aggregated" src="https://user-images.githubusercontent.com/8429802/51890650-dcee8400-239c-11e9-8b62-9f9009950e5e.png">
