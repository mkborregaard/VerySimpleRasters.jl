module VerySimpleRasters

using HDF5
using RecipesBase
using Mmap
using Dates

import Base.getindex

include("datatypes.jl")
include("importASCII.jl")
include("grd_file.jl")
include("operations.jl")
include("plotrecipes.jl")

export VerySimpleRaster, importASCII, write_grd
export extract, aggregate, crop, mask
export writeraster
end # module
