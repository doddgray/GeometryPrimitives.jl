# Assume that the prism axis is orthogonal to the prism base.  This makes it much easier to
# implement surfpt_nearby by separating the axis and base dimensions (i.e., the closest
# point on the surface can be found independently in axis and base dimensions).
#
# Supporting skewed prisms is not impossible, but this will be done in the future if the
# demand is high.

export Prism

struct Prism{B<:Shape2,T<:Number} <: Shape3
    c::SVector{3,T}  # prism center
    b::B  # base shape described in prism coordinates (i.e, when translating prism, do not need to translate b)
    h2::T  # height * 0.5
    p::SMatrix{3,3,T,9}  # projection matrix to prism coordinates; must be orthonormal (see surfpt_nearby)
    Prism{B,T}(c,b,h2,p) where {B,T} = new(c,b,h2,p)  # suppress default outer constructor
end

Prism{B}(c::SVector{3,<:Number}, b::B, h2::Number, p::SMatrix{3,3,<:Number}) where {B<:Shape2} =
    (T = promote_eltype(eltype(c), typeof(h2), eltype(p)); Prism{B,T}(c, b, h2, p))

Prism(c::SVector{3,<:Number},
      b::B,
      h::Number=Inf,
      axes::SMatrix{3,3,<:Number,9}=SMatrix{3,3,Float64}(I)  # columns are axes vectors: first two columns span prism base, and last column is prism axis
      ) where {B<:Shape2} =
    Prism{B}(c, b, 0.5h, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))))

Prism(c::AbstractVector{<:Number}, b::Shape2, h::Number=Inf, axes::AbstractMatrix{<:Number}=Matrix{Float64}(I,length(c),length(c))) =
    Prism(SVector{3}(c), b, h, SMatrix{3,3}(axes))

Base.:(==)(s1::Prism, s2::Prism) = s1.c==s2.c && s1.b==s2.b && s1.h2==s2.h2 && s1.p==s2.p
Base.isapprox(s1::Prism, s2::Prism) = s1.c≈s2.c && s1.b≈s2.b && s1.h2≈s2.h2 && s1.p≈s2.p
Base.hash(s::Prism, h::UInt) = hash(s.c, hash(s.b, hash(s.h2, hash(s.p, hash(:Prism, h)))))

translate(s::Prism{B}, ∆::SVector{3,<:Number}) where {B} = Prism{B}(s.c + ∆, s.b, s.h2, s.p)

function level(x::SVector{3,<:Number}, s::Prism)
    y = s.p * (x - s.c)  # coordinates after projection
    ya = y[3]  # scalar: coordinate in axis dimension
    yb = y[SVector(1,2)]  # SVector{2}: coordinate in base dimensions

    return min(1 - abs(ya)/s.h2, level(yb,s.b))
end

function surfpt_nearby(x::SVector{3,<:Number}, s::Prism)
    ax = s.p'  # prism axes: columns are not only unit vectors, but also orthogonal

    y = s.p * (x - s.c)  # x in prism coordinates
    ya = y[3]  # scalar: coordinate in axis dimension
    yb = y[SVector(1,2)]  # SVector{2}: coordinates in base dimensions

    la = abs(ya)
    abs∆a = abs(s.h2 - la)  # scalar: distance between x and base point closest to x
    surfa = SVector(yb[1], yb[2], copysign(s.h2, ya))  # SVector{3}: coordinates of base point closest to x
    nouta = SVector(zero(ya), zero(ya), copysign(one(ya), ya))  # SVector{3}: outward direction normal at surfa
    onbnda = abs∆a ≤ rtol(s.h2)
    isouta = s.h2<la || onbnda

    surfb2, noutb2 = surfpt_nearby(yb, s.b)  # (SVector{2}, SVector{2}): side point closest to x and outward direction normal to side there
    abs∆b = norm(surfb2 - yb)  # scalar: distance between x and side point closest to x
    surfb = SVector(surfb2[1], surfb2[2], ya)  # SVector{3}: coordinates of side point closest to x
    noutb = SVector(noutb2[1], noutb2[2], zero(eltype(noutb2)))  # SVector{3}: outward direction normal to side surface at surfb
    basesize = abs.((-)(bounds(s.b)...))  # SVector{2}: size of bounding rectancle of base
    onbndb = abs∆b ≤ rtol(maximum(basesize))
    isoutb = yb∉s.b || onbndb

    if isouta && isoutb  # x outside in both axis and base dimensions
        surf = SVector(surfb[1], surfb[2], surfa[3])
        nout = (onbnda && onbndb) ? (noutb + nouta) : (y - surf)
        nout = norm(nout)==Inf ? isinf.(nout) .* sign.(nout) : normalize(nout)  # e.g., return [0,0,-1] for nout = [1,-2,-Inf]
    elseif !isouta && isoutb  # x outside in base dimensions, but inside prism in axis dimension
        (surf, nout) = (surfb, noutb)
    elseif isouta && !isoutb # x outside in axis dimension, but inside prism in base dimensions
        (surf, nout) = (surfa, nouta)
    else  # !isouta && !isoutb: x strictly inside prism
        (surf, nout) = (abs∆a ≤ abs∆b) ? (surfa, nouta) : (surfb, noutb)
    end

    # Transform back from prism coordinates: x = ax*y + s.c.  (The upstream expression
    # ax*(surf + s.c) is equivalent only when ax = I or s.c = 0.)
    return ax*surf + s.c, ax*nout
end

function bounds(s::Prism)
    ax = s.p'  # prism axes: columns are not only unit vectors, but also orthogonal
    a = ax[:,3]  # SVector{3}
    h2a = s.h2 * a

    l0, u0 = bounds_ctrcut(s)  # (SVector{3}, SVector{3})
    l1, u1 = l0+h2a, u0+h2a  # (SVector{3}, SVector{3})
    l2, u2 = l0-h2a, u0-h2a  # (SVector{3}, SVector{3})

    return min.(l1,l2)+s.c, max.(u1,u2)+s.c
end


include("cylinder.jl")
include("polygonal.jl")
include("sectoral.jl")
