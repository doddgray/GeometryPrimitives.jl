export Sphere

mutable struct Sphere{N,N²,D,T} <: Shape{N,N²,D,T}
    c::SVector{N,T}  # center of sphere
    r::T  # radius
    data::D  # auxiliary data
    Sphere{N,N²,D,T}(c,r,data) where {N,N²,D,T<:Real} = new(c,r,data)  # suppress default outer constructor
end

Sphere(c::SVector{N,T}, r::T, data::D=nothing) where {N,D,T<:Real} = Sphere{N,N*N,D,T}(c, r, data)
Sphere(c::AbstractVector{<:Real}, r::<:Real, data=nothing) = (N = length(c); Sphere(SVector{N}(c), r, data))

Base.:(==)(s1::Sphere, s2::Sphere) = s1.c==s2.c && s1.r==s2.r && s1.data==s2.data
Base.isapprox(s1::Sphere, s2::Sphere) = s1.c≈s2.c && s1.r≈s2.r && s1.data==s2.data
Base.hash(s::Sphere, h::UInt) = hash(s.c, hash(s.r, hash(s.data, hash(:Sphere, h))))

Base.in(x::SVector{N,<:Real}, s::Sphere{N}) where {N} = sum(abs2, x - s.c) ≤ s.r^2

function surfpt_nearby(x::SVector{N,<:Real}, s::Sphere{N}) where {N}
    nout = x==s.c ? SVector(ntuple(k -> k==1 ? 1.0 : 0.0, Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

translate(s::Sphere{N,N²,D}, ∆::SVector{N,T}) where {N,N²,D,T<:Real} = Sphere{N,N²,D,T}(s.c+∆, s.r, s.data)

bounds(s::Sphere) = (s.c.-s.r, s.c.+s.r)
