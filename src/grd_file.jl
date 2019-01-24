function read_grd(filename)
    a = open(filename, "r") do IO
        read(IO, String)
    end
    b = filter!(x->!isempty(x) && !(x[1] == '['), split(a, r"\r\n"))
    Dict(x => replace(y, " " => "") for (x,y) in split.(b, Ref("=")))
end

#data = read_grd("data/Gomera_dem.grd")
#using Mmap

#b = open("data/Gomera_dem.gri", "r")
#mm = Mmap.mmap(b, Matrix{Float32}, (167, 127))
#close(b)


function VerySimpleRaster(grdfile::String)
    @assert grdfile[end-3:end] ∈ (".grd", ".gri")
    fname = grdfile[1:end-4]
    data = read_grd(fname*".grd")
    ncols = parse(Int, data["nrows"])
    nrows = parse(Int, data["ncols"])
    dt = datatype_translation[data["datatype"]]
    fname = abspath(fname*".gri")
    mat = open(fname, "r") do IO
        Mmap.mmap(IO, Matrix{dt}, (nrows, ncols))
    end
    xmin, xmax = parse(Float64, data["xmin"]), parse(Float64, data["xmax"])
    ymin, ymax = parse(Float64, data["ymin"]), parse(Float64, data["ymax"])
    celly = (ymax - ymin) / ncols
    cellx = (xmax - xmin) / nrows
    VerySimpleRaster(mat, parse(dt, data["nodatavalue"]), xmin, ymin, cellx, celly, data["projection"], fname)
end

datatype_translation = Dict{String, DataType}(
    "LOG1S" => Bool,
    "INT1S" => Int8,
    "INT2S" => Int16,
    "INT4S" => Int32,
    "INT8S" => Int64,
    "INT1U" => UInt8,
    "INT2U" => UInt16,
    "FLT4S" => Float32,
    "FLT8S" => Float64
)
