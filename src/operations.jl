function aggregate(mymatrix, factor::Int, fun::Function)
    newm, newn = ceil.(Int,size(mymatrix)./factor)
    ret = Matrix{typeof(fun([first(mymatrix)]))}(undef, newm, newn)
    for j in 1:newn-1
       for i in 1:newm-1
          ret[i,j] = fun(vec(mymatrix[(i-1)*factor .+ (1:factor), (j-1)*factor .+ (1:factor)]))
       end
       ret[newm,j] = fun(vec(mymatrix[(newm-1)*factor + 1:size(mymatrix,1), (j-1)*factor .+ (1:factor)]))
    end
    for i in 1:newm-1
       ret[i, newn] = fun(vec(mymatrix[(i-1)*factor .+ (1:factor), (newn-1)*factor + 1:size(mymatrix,2)]))
    end
    ret[newm, newn] = fun(vec(mymatrix[(newm-1)*factor + 1:size(mymatrix, 1), (newn-1)*factor + 1:size(mymatrix,2)]))
    ret
end

extract(vsr::VerySimpleRaster, x, y) = extract(vsr, (x, y))
function extract(vsr::VerySimpleRaster, tup)
   x, y = tup[1], tup[2]
   indx, indy = findfirst(i-> i>x, vsr.xs), findfirst(i-> i>y, vsr.ys)
   vsr[indx, indy]
end

extract(vsr::VerySimpleRaster, tup::Tuple{T, T}) where T <: AbstractVector = extract.(Ref(vsr), tup)
extract(vsr::VerySimpleRaster, ar::AbstractVector) = extract.(Ref(vsr), ar)


function getindex(vsr::VerySimpleRaster, x...)
   ret = vsr.mat[x...]
   ret == vsr.nodata && return missing
   ret
end

function getindex(vsr::VerySimpleRaster, x)
   ret = vsr.mat[x]
   ret == vsr.nodata && return missing
   ret
end
