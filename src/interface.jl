eachcoordinate(vsr::VerySimpleRaster) = product(vsr.xs, vsr.ys)
coordinates(vsr::VerySimpleRaster) = collect(flatten(eachcoordinate(vsr)))
