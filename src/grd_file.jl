function read_grd(filename)
    a = readlines(filename)
    filter!(x->!isempty(x) && !(x[1] == '['), a)
    Dict((c = match(r"([^=]+)=(.*)", st); c.captures[1] => strip(c.captures[2])) for st in a)
end

function write_grd(grdfile::String, vsr::VerySimpleRaster)
    grdfile[end-3:end] ∈ (".grd", ".gri") && (grdfile = grdfile[1:end-4])
    open(grdfile*".gri", "w") do IO
        write(IO, vsr.mat)
    end
    write_grd_header(grdfile, vsr)
end

function write_grd_header(grdfile, vsr::VerySimpleRaster)
    bb = bbox(vsr)
    write_grd_header(grdfile,
        nrows= size(vsr.mat, 2),
        ncols= size(vsr.mat, 1),
        xmin= bb.xmin,
        ymin= bb.ymin,
        xmax= bb.xmax,
        ymax= bb.ymax,
        projection= vsr.projection,
        datatype= eltype(vsr.mat),
        nodata= vsr.nodata,
        minvalue = minimum(filter(x->x!=vsr.nodata, vsr.mat)),
        maxvalue= maximum(filter(x->x!=vsr.nodata, vsr.mat))
    )
end

function write_grd_header(grdfile;
                        nrows = error("unspecified argument nrows"), ncols = error("unspecified argument ncols"),
                        xmin = error("unspecified argument xmin"), ymin = error("unspecified argument ymin"),
                        xmax = error("unspecified argument xmax"), ymax = error("unspecified argument ymax"),
                        projection = error("unspecified argument projection"), datatype = error("unspecified argument datatype"),
                        nodata = error("unspecified argument nodata"),
                        minvalue = error("unspecified argument minvalue"), maxvalue = error("unspecified argument maxvalue"))
    grdfile[end-3:end] == ".grd" || (grdfile = grdfile * ".grd")
    a = open(grdfile, "w") do IO
        write(IO, """
        [general]
        creator=Julia Package VerySimpleRasters
        created= $(string(now()))
        [georeference]
        nrows= $(nrows)
        ncols= $(ncols)
        xmin= $(xmin)
        ymin= $(ymin)
        xmax= $(xmax)
        ymax= $(ymax)
        projection= $(projection)
        [data]
        datatype= $(rev_datatype_translation[datatype])
        nodatavalue= $(nodata)
        byteorder= little
        nbands= 1
        minvalue= $(minvalue)
        maxvalue= $(maxvalue)
        [description]
        layername= layer1""")
    end
end


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

const datatype_translation = Dict{String, DataType}(
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

const rev_datatype_translation = Dict{DataType, String}(v => k for (k,v) in datatype_translation)

"""
    writeraster(file::String, vsr::VerySimpleRaster)

Writes the raster as an R .grd file. Under the hood simply copies the temporary
on-disk version of the raster to `file`.
"""
function writeraster(fname::String, vsr::VerySimpleRaster)
    fname[end-3:end] ∈ (".gri", "grd") && (fname = fname[1:end-4])
    cp(vsr.filename, fname*".gri")
    cp(vsr.filename[end-3:end]*".grd", fname*".grd")
end
