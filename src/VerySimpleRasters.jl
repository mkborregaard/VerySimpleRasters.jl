module VerySimpleRasters

using HDF5
using RecipesBase
using Mmap

include("datatypes.jl")
include("importASCII.jl")
include("plotrecipes.jl")

export VerySimpleRaster, importASCII
end # module
