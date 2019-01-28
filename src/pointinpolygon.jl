
 struct Pt
     x::Float64
     y::Float64
end

function isinside(r::Pt, poly::AbstractVector{<:Pt})
    # An implementation of Hormann-Agathos (2001) Point in Polygon algorithm
    # See: http://www.sciencedirect.com/science/article/pii/S0925772101000128
    # Code segment adapted from PolygonClipping.jl
    c = false
    detq(q1,q2,r) = (q1[1]-r.x)*(q2[2]-r.y)-(q2[1]-r.x)*(q1[2]-r.y)
    for i in eachindex(poly)[2:end]
        q2 = poly[i]
        q1 = poly[i-1]
        if q1 == r
            @warn("point on polygon vertex - returning false")
            return false
        end
        if q2.y == r.y
            if q2.x == r.x
                @warn("point on polygon vertex - returning false")
                return false
            elseif (q1.y == r.y) && ((q2.x > x) == (q1.x < r.x))
                @warn("point on edge - returning false")
                return false
            end
        end
        if (q1.y < r.y) != (q2.y < r.y) # crossing
            if q1.x >= r.x
                if q2.x > r.x
                    c = !c
                elseif ((detq(q1,q2,r) > 0) == (q2.y > q1.y)) # right crossing
                    c = !c
                end
            elseif q2.x > r.x
                if ((detq(q1,q2,r) > 0) == (q2.y > q1.y)) # right crossing
                    c = !c
                end
            end
        end
    end
    return c
end













# # fast algorithm from https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
#
# struct Pt
#     x::Float64
#     y::Float64
# end
#
# # takes two points and generates a function giving the line. All points giving 0 from that function will be on the line
# _points_to_line(pt1::Pt, pt2::Pt) = (pt::Pt) -> (pt2.y - pt1.y)*pt.x + (pt2.x - pt1.x)*pt.y + (pt2.x*pt1.y - pt1.x*pt2.y)
# function _lineintersects(pt1_1, pt1_2, pt2_1, pt2_2,
#                         linefun1 = _points_to_line(pt1_1, pt1_2), linefun2 = _points_to_line(pt2_1, pt2_2))
#     (pt1_1 == pt1_2 || pt2_1 == pt2_2) && return 0 # line has length 0
#     tmp1, tmp2 = linefun1(pt2_1), linefun1(pt2_2)
#     sign(tmp1) == sign(tmp2) && return 0
#     tmp3, tmp4 = linefun2(pt1_1), linefun2(pt1_2)
#     sign(tmp3) == sign(tmp4) && return 0
#     return 1
# end
#
#
# function point_in_polygon(pt::Pt, polygon::AbstractVector{<:Pt},
#                         bbox = (extrema(p.x for p in polygon)..., extrema(p.y for p in polygon)...))
#
#     (pt.x < bbox[1] || pt.x > bbox[2] || pt.y < bbox[3] || pt.y > bbox[4] || bbox[1] == bbox[2] || bbox[3] == bbox[4]) && return false
#
#     raypt = Pt(bbox[1] - (bbox[2] - bbox[1]), pt.y)
#     rayline = _points_to_line(raypt, pt)
#     intersections = 0
#
#     @inbounds for i in eachindex(polygon)[2:end]
#         intersections += _lineintersects(raypt, pt, polygon[i-1], polygon[i], rayline)
#     end
#     intersections += _lineintersects(raypt, pt, last(polygon), first(polygon), rayline)
#     Bool(intersections%2)
# end
