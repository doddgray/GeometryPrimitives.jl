export Sphere

mutable struct Sphere{N,N²,D,T} <: Shape{N,N²,D,T}
    c::SVector{N,T}  # center of sphere
    r::T  # radius
    data::D  # auxiliary data
    # Sphere{N,N²,D,T}(c,r,data) where {N,N²,D,T<:Real} = new(c,r,data)  # suppress default outer constructor
end

Sphere(c::SVector{N,T}, r::T, data::D=nothing) where {N,D,T<:Real} = Sphere{N,N*N,D,T}(c, r, data)
Sphere(c::AbstractVector{<:Real}, r::Real, data=nothing) = (N = length(c); Sphere(SVector{N}(c), r, data))

Base.:(==)(s1::Sphere, s2::Sphere) = s1.c==s2.c && s1.r==s2.r && s1.data==s2.data
Base.isapprox(s1::Sphere, s2::Sphere) = s1.c≈s2.c && s1.r≈s2.r && s1.data==s2.data
Base.hash(s::Sphere, h::UInt) = hash(s.c, hash(s.r, hash(s.data, hash(:Sphere, h))))

Base.in(x::SVector{N,<:Real}, s::Sphere{N}) where {N} = sum(abs2, x - s.c) ≤ s.r^2

function surfpt_nearby(x::SVector{N,<:Real}, s::Sphere{N}) where {N}
    nout = x==s.c ? SVector(ntuple(k -> k==1 ? 1.0 : 0.0, Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

translate(s::Sphere{N,N²,D}, ∆::SVector{N}) where {N,N²,D,T<:Real} = Sphere{N,N²,D,T}(s.c+∆, s.r, s.data)

bounds(s::Sphere) = (s.c.-s.r, s.c.+s.r)


"""
ChainRulesCore differentiation rules for Sphere

Rules for reverse-mode differentiation (`rrules`) below just let the AD system auto-differentiate the Sphere
contstructor and `surfpt_nearby` method for Spheres, but then re-format the gradient result as a `Tangent{Sphere}`,
Sphere-like struct with addition methods defined, enabling gradient accumulation.
"""

function rrule(TS::Type{<:Sphere}, c, r, data)
    Sphere_pullback(Δsphere) = NoTangent(), Δsphere.c, Δsphere.r, Δsphere.data
    return TS(c,r,data), Sphere_pullback
end

function rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(surfpt_nearby), x::SVector{N,T}, s::Sphere{N}) where {N,T<:Real}
    r₀_n⃗_and_surfpt_nearby_Sphere_fields_pullback = rrule_via_ad(config,(x,c,r)->surfpt_nearby(x,Sphere(c,r,nothing)),x,s.c,s.r)
    r₀_n⃗    = first(r₀_n⃗_and_surfpt_nearby_Sphere_fields_pullback)
    surfpt_nearby_Sphere_fields_pullback = last( r₀_n⃗_and_surfpt_nearby_Sphere_fields_pullback )
    function surfpt_nearby_Sphere_pullback(r₀_bar,n⃗_bar)
        x̄,c̄,r̄ = surfpt_nearby_Sphere_fields_pullback(r₀_bar,n⃗_bar)
        return NoTangent(), x̄, @thunk(canonicalize(Tangent{typeof(s)}(;c=c̄,r=r̄)))
    end
    return r₀_n⃗, surfpt_nearby_Sphere_pullback
end