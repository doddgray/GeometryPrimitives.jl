export Ellipsoid

struct Ellipsoid{N,N²,T<:Number} <: Shape{N,N²}
    c::SVector{N,T}  # center of ellipsoid
    ri2::SVector{N,T}  # inverse squares of "radii" (semi-axes) in axis directions
    p::SMatrix{N,N,T,N²}  # projection matrix to Ellipsoid coordinates; must be orthonormal (see surfpt_nearby)
    Ellipsoid{N,N²,T}(c,ri2,p) where {N,N²,T} = new(c,ri2,p)  # suppress default outer constructor
end

Ellipsoid{N,N²}(c::SVector{N,<:Number}, ri2::SVector{N,<:Number}, p::SMatrix{N,N,<:Number}) where {N,N²} =
    (T = promote_eltype(eltype(c), eltype(ri2), eltype(p)); Ellipsoid{N,N²,T}(c, ri2, p))

Ellipsoid(c::SVector{N,<:Number},
          r::SVector{N,<:Number},
          axes::SMatrix{N,N,<:Number}=SMatrix{N,N,Float64}(I)
          ) where {N} =
    Ellipsoid{N,N*N}(c, inv.(r).^2, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))))

Ellipsoid(c::AbstractVector{<:Number},  # center of ellipsoid
          r::AbstractVector{<:Number},  # ""
          axes::AbstractMatrix{<:Number}=Matrix{Float64}(I,length(c),length(c))) =  # columns are axes vector; assumed orthogonal
    (N = length(c); Ellipsoid(SVector{N}(c), SVector{N}(r), SMatrix{N,N}(axes)))

Ellipsoid(s::Cuboid{N,N²}) where {N,N²} = Ellipsoid{N,N²}(s.c, (s.r).^-2, s.p)

Base.:(==)(s1::Ellipsoid, s2::Ellipsoid) = s1.c==s2.c && s1.ri2==s2.ri2 && s1.p==s2.p
Base.isapprox(s1::Ellipsoid, s2::Ellipsoid) = s1.c≈s2.c && s1.ri2≈s2.ri2 && s1.p≈s2.p
Base.hash(s::Ellipsoid, h::UInt) = hash(s.c, hash(s.ri2, hash(s.p, hash(:Ellipsoid, h))))

translate(s::Ellipsoid{N,N²}, ∆::SVector{N,<:Number}) where {N,N²} = Ellipsoid{N,N²}(s.c + ∆, s.ri2, s.p)

# Written with broadcasting and sum() instead of dot() because traced array types
# (Reactant) support elementwise broadcasting of static arrays but not whole-array
# arithmetic or dot().
level(x::SVector{N,<:Number}, s::Ellipsoid{N}) where {N} = 1 - √(sum((s.p * (x .- s.c)).^2 .* s.ri2))

function surfpt_nearby(x::SVector{N,<:Number}, s::Ellipsoid{N}) where {N}
    if x == s.c
        _m, i = findmax(s.ri2)
        nout = s.p[i,:]  # assume s.p is orthogonal
        return s.c + nout/√s.ri2[i], nout
    end

    # For a given point x and equation of ellipsoid f(x) = 1, find t such that x₀ = x + t*∇f(x)
    # is on the ellipsoid.  Eventually this reduces to a quadratic equation for t.  The
    # following is evaluation of the quadratic formula in a numerically stable way.
    px = s.p * (x - s.c)  # in ellipsoid's coordinates
    px2 = px.^2

    px²r⁻² = px2 ⋅ s.ri2
    px²r⁻⁴ = px2 ⋅ (s.ri2.^2)
    px²r⁻⁶ = px2 ⋅ (s.ri2.^3)

    q24 = (px²r⁻² - 1) / px²r⁻⁴
    q64 = px²r⁻⁶ / px²r⁻⁴

    t = -q24 / (1 + √(1 - q24 * q64))

    # From t, recover x₀ = x + t*∇f(x).
    px₀ = (t*s.ri2 .+ 1) .* px  # surface point in ellipsoid coordinates

    # Transform back to the original coordinates.
    x₀ = s.p' * px₀ + s.c
    nout = normalize(s.p' * (px .* s.ri2))

    return x₀, nout
end

function boundpts(s::Ellipsoid{N}) where {N}
    # Return the points tangential to the bounding box.
    # For N = 3, it returns three points at which the direction normals are +x, +y, +z
    # directions, respectively.

    r2 = 1 ./ s.ri2
    ndir = s.p  # Cartesian directions in ellipsoid coordinates: s.p * I

    # In the ellipsoid coordinates, the point on the ellipsoid where the direction normal is
    # n is (n .* r2) / sqrt(r2' * n.^2).  Below is the broadcasted version of this over a
    # matrix n, whose each column is a direction normal.  Once calculated, we need to
    # change the coordinates back to the original coordinates.
    M = s.p' * ((ndir .* r2) ./ sqrt.(r2' * ndir.^2))

    return M
end

function bounds(s::Ellipsoid{N}) where {N}
    M = boundpts(s)

    # Note that when M does not include NaN, we can simply set m = diag(M), because the
    # first (second, third) column of M is the point on the ellipsoid at which the normal
    # vector is +x (+y, +z).  Therefore, the x (y, z) point of the first (second, third)
    # column has the largest x (y, z) coordinate.
    #
    # However, if one of a, b, c is zero, the shape is a disk.  Then one column of M can be
    # completely filled with NaN.  This column must not be counted as a bounding point, so
    # we apply NaN-ignoring maximum by StaticArrays.reducedim along the row direction.
    #
    # For the efficient implementation of NaN-ignoring maximum, see
    # https://discourse.julialang.org/t/inconsistency-between-findmax-and-maximum-with-respect-to-nan/4214/8
    m = reduce((x,y) -> x<y ? y : x, M, init=-Inf, dims=Val(2))[:,1]

    return (s.c-m,s.c+m)
end
