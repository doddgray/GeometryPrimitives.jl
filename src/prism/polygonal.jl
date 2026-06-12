export PolygonalPrism

const PolygonalPrism{K,K2,T} = Prism{Polygon{K,K2,T},T}

# Below, if we called PolygonalPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Polygon{K,K2,T},T}(c, ...), which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of PolygonalPrism(c, ...).
function PolygonalPrism(c::SVector{3,<:Real},
                        v::SMatrix{2,K,<:Real},  # 2D coordinates of base vertices in projected prism coordinates
                        h::Real=Inf,
                        a::SVector{3,<:Real}=SVector(0.0,0.0,1.0)
                        ) where {K}
    T = promote_eltype(eltype(c), eltype(v), typeof(h), eltype(a))
    â = normalize(SVector{3,T}(a))
    return Prism(SVector{3,T}(c), Polygon(SMatrix{2,K,T}(v)), h, [orthoaxes(â)... â])
end

PolygonalPrism(c::AbstractVector{<:Real},  # center of prism
               v::AbstractMatrix{<:Real},  # vertices of base polygon
               h::Real=Inf,  # height of prism
               a::AbstractVector{<:Real}=[0.0,0.0,1.0]) =  # axis direction of prism
    (K = size(v,1); PolygonalPrism(SVector{3}(c), SMatrix{2,K}(v), h, SVector{3}(a)))

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::Prism{<:Polygon{K}}) where {K}
    p = s.p'  # projection matrix to prism coordinates: rows are not only unit vectors, but also orthogonal
    v = [s.b.v; @SMatrix(zeros(1,K))]  # SMatrix{3,K}: 3D vectices in prism axis coordinates
    w = p * v  # SMatrix{3,K}: vertices in external coordinates

    return minimum(w, dims=Val(2))[:,1], maximum(w, dims=Val(2))[:,1]
end
