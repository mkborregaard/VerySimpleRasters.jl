module VerySimpleRasters

using HDF5
using RecipesBase

include("datatypes.jl")
include("importASCII.jl")
include("plotrecipes.jl")

export VerySimpleRaster, importASCII
end # module
