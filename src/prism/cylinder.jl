export Cylinder

const Cylinder = Prism{Sphere{2,4,Nothing}}

# Below, if we called Cylinder(c, ...) in the function body, it would call the inner
# constructor Prism{Sphere{2,4,Nothing}}(c, ...) because Cylinder = Prism{Sphere{2,4,Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of Cylinder(c, ...).
function Cylinder(c::SVector{3,T}, r::T, h::T=typemax(T), a::SVector{3,T}=SVector{3}(0.0,0.0,1.0), data::D=nothing) where {T<:Real,D}
    â = normalize(a)
    return Prism(c, Sphere(SVector{2,T}(0.0,0.0),r, nothing), h, [orthoaxes(â)... â], data)
end

# Cylinder(c::AbstractVector{<:Real},  # center of cylinder
#          r::Real,  # radius of base
#          h::Real=Inf,  # height of cylinder
#          a::AbstractVector{<:Real}=[0.0,0.0,1.0],  # axis direction of cylinder
#          data=nothing) =
#     Cylinder(SVector{3}(c), r, h, SVector{3}(a), data)

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::Cylinder)
    ax = inv(s.p)  # prism axes: columns are not only unit vectors, but also orthogonal
    r = s.b.r
    el = Ellipsoid(SVector{3}(0.0,0.0,0.0), SVector{3}(r,r,0.0), ax)  # center is set at origin to return bounds with respect to prism center

    return bounds(el)
end
