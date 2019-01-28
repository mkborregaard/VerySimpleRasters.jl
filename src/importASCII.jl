
function importASCII(x::AbstractString, outfile = "")
    file = open(x, "r")
    nc = parse(Int, match(r"NCOLS (.+)", readline(file)).captures[1])
    nr = parse(Int, match(r"NROWS (.+)", readline(file)).captures[1])
    xll = parse(Float64, match(r"XLLCORNER (.+)", readline(file)).captures[1])
    yll = parse(Float64, match(r"YLLCORNER (.+)", readline(file)).captures[1])
    cell = parse(Float64, match(r"CELLSIZE (.+)", readline(file)).captures[1])
    NA = parse(Float64, match(r"NODATA_value (.+)", readline(file)).captures[1])
    tmp = Vector{Float64}(undef, nc)

    outfile == "" && (outfile = tempname())
    outfile[end-3:end] ∈ (".gri", "grd") && (outfile = outfile[1:end-4])
    IO = open(outfile*".gri", "w")

    for row in nr:-1:1
        tmp .= parse.(Float64, split(readline(file), " "))
        for r in tmp
            write(IO, r)
        end
    end
    close(file)
    close(IO)
    mat = Mmap.mmap(outfile*".gri", Matrix{Float64}, (nc, nr))
    vsr = VerySimpleRaster(mat, NA, xll, yll, cell, cell, "", outfile*".gri")
    write_grd_header(outfile*".grd", vsr)
    vsr
end
