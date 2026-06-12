export Ball

struct Ball{N,N²,T<:Number} <: Shape{N,N²}
    c::SVector{N,T}  # center of ball
    r::T  # radius
end

Ball{N,N²}(c::SVector{N,<:Number}, r::Number) where {N,N²} =
    (T = promote_eltype(eltype(c), typeof(r)); Ball{N,N²,T}(c, r))
Ball(c::SVector{N,<:Number}, r::Number) where {N} = Ball{N,N*N}(c, r)
Ball(c::AbstractVector{<:Number}, r::Number) = (N = length(c); Ball(SVector{N}(c), r))

Base.:(==)(s1::Ball, s2::Ball) = s1.c==s2.c && s1.r==s2.r
Base.isapprox(s1::Ball, s2::Ball) = s1.c≈s2.c && s1.r≈s2.r
Base.hash(s::Ball, h::UInt) = hash(s.c, hash(s.r, hash(:Ball, h)))

translate(s::Ball{N}, ∆::SVector{N,<:Number}) where {N} = Ball(s.c + ∆, s.r)

# Written with broadcasting (x .- s.c) because traced array types (Reactant) support
# elementwise broadcasting of static arrays but not whole-array arithmetic.
level(x::SVector{N,<:Number}, s::Ball{N}) where {N} = 1 - √(sum(abs2, x .- s.c) / s.r^2)

function surfpt_nearby(x::SVector{N,<:Number}, s::Ball{N,N²,T}) where {N,N²,T}
    nout = x==s.c ? SVector(ntuple(k -> k==1 ? one(T) : zero(T), Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

bounds(s::Ball) = (s.c.-s.r, s.c.+s.r)
