export Cylinder

const Cylinder{T} = Prism{Ball{2,4,T},T}

# Below, if we called Cylinder(c, ...) in the function body, it would call the inner
# constructor Prism{Ball{2,4,T},T}(c, ...), which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of Cylinder(c, ...).
function Cylinder(c::SVector{3,<:Number},
                  r::Number,
                  h::Number=Inf,
                  a::SVector{3,<:Number}=SVector(0.0,0.0,1.0))
    T = promote_eltype(eltype(c), typeof(r), typeof(h), eltype(a))
    â = normalize(SVector{3,T}(a))
    return Prism(SVector{3,T}(c), Ball(SVector(zero(T),zero(T)), r), h, [orthoaxes(â)... â])
end

Cylinder(c::AbstractVector{<:Number},  # center of cylinder
         r::Number,  # radius of base
         h::Number=Inf,  # height of cylinder
         a::AbstractVector{<:Number}=[0.0,0.0,1.0]) =  # axis direction of cylinder
    Cylinder(SVector{3}(c), r, h, SVector{3}(a))

# Return the bounds of the center cut with respect to the prism center.
# (Dispatch on the base-shape type rather than the Cylinder alias so that prisms whose
# base and prism element types differ are also covered.)
function bounds_ctrcut(s::Prism{<:Ball{2}})
    ax = s.p'  # prism axes: columns are not only unit vectors, but also orthogonal
    r = s.b.r
    el = Ellipsoid(SVector(zero(r),zero(r),zero(r)), SVector(r,r,zero(r)), ax)  # center is set at origin to return bounds with respect to prism center

    return bounds(el)
end
