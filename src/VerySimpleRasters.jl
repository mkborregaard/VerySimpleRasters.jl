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

export VerySimpleRaster, importASCII
export extract, aggregate
export writeraster
end # module
