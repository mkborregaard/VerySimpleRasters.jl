
function importASCII(x::AbstractString)
    file = open(x, "r")
    nc = parse(Int, match(r"NCOLS (.+)", readline(file)).captures[1])
    nr = parse(Int, match(r"NROWS (.+)", readline(file)).captures[1])
    xll = parse(Float64, match(r"XLLCORNER (.+)", readline(file)).captures[1])
    yll = parse(Float64, match(r"YLLCORNER (.+)", readline(file)).captures[1])
    cell = parse(Float64, match(r"CELLSIZE (.+)", readline(file)).captures[1])
    NA = parse(Float64, match(r"NODATA_value (.+)", readline(file)).captures[1])
    ret = Matrix{Float64}(undef, nr, nc)
    for row in nr:-1:1
        ret[row,:] .= parse.(Float64, split(readline(file), " "))
    end
    close(file)
    VerySimpleRaster(ret, xll, yll, cell, NA)
end
