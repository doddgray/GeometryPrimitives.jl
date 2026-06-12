export Polygon

#= Polygon (for a base shape) =#

# Assume the followings for the polygon represented by Polygon:
# - The polygon is convex.
# - The vertices are listed in the counter-clockwise order around the origin.
struct Polygon{K,K2,T<:Number} <: Shape2  # K2 = 2K
    v::SMatrix{2,K,T,K2}  # vertices
    n::SMatrix{2,K,T,K2}  # direction normals to edges
    Polygon{K,K2,T}(v,n) where {K,K2,T} = new(v,n)  # suppress default outer constructor
end

Polygon{K,K2}(v::SMatrix{2,K,<:Number}, n::SMatrix{2,K,<:Number}) where {K,K2} =
    (T = promote_eltype(eltype(v), eltype(n)); Polygon{K,K2,T}(v, n))

function Polygon(v::SMatrix{2,K,<:Number}) where {K}
    # Sort the vertices in the counter-clockwise direction
    w = v .- mean(v, dims=Val(2))  # v in center-of-mass coordinates
    ϕ = mod.(atan.(w[2,:], w[1,:]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
    if !issorted(ϕ)
        # Do this only when ϕ is not sorted, because the following uses allocations.
        ind = SVector{K}(sortperm(ϕ))
        v = v[:,ind]  # SMatrix{2,K}: sorted v
    end

    # Calculate the increases in angle between neighboring edges.
    ∆v = hcat(diff(v, dims=Val(2)), SMatrix{2,1}(v[:,1]-v[:,end]))  # SMatrix{2,K}: edge directions
    icurr = ntuple(identity, Val(K-1))
    inext = ntuple(x->x+1, Val(K-1))

    # ∆ϕ[k] is the angle from the edge direction ∆v[:,k] to ∆v[:,k+1], calculated as
    # atan(∆v[:,k] × ∆v[:,k+1], ∆v[:,k] ⋅ ∆v[:,k+1]) ∈ (-π, π].  (This is the angle of the
    # complex number ratio ∆z[k+1] / ∆z[k] for ∆z = ∆v[1,:] + im*∆v[2,:], written without
    # complex arithmetic to remain differentiable by the widest range of AD tools.)
    ∆vc = ∆v[:,SVector(icurr)]  # SMatrix{2,K-1}: current edge directions
    ∆vn = ∆v[:,SVector(inext)]  # SMatrix{2,K-1}: next edge directions
    ∆ϕ = atan.(∆vc[1,:].*∆vn[2,:] .- ∆vc[2,:].*∆vn[1,:], ∆vc[1,:].*∆vn[1,:] .+ ∆vc[2,:].*∆vn[2,:])

    # Check all the angle increases are positive.  If they aren't, the polygon is not convex.
    all(∆ϕ .> 0) || throw("v = $v should represent vertices of convex polygon.")

    n = [∆v[2,:] -∆v[1,:]]'  # SMatrix{2,K}; outward normal directions to edges
    n = n ./ hypot.(n[1,:], n[2,:])'  # normalize

    return Polygon{K,2K}(v,n)
end

Polygon(v::AbstractMatrix{<:Number}) = (K = size(v,2); Polygon(SMatrix{2,K}(v)))

# Regular polygon
function Polygon{K}(c::SVector{2,<:Number},
                    r::Number,  # distance between center and each vertex
                    θ::Number=0.0  # angle from +y-direction towards first vertex
                    ) where {K}
    ∆θ = 2π / K

    θs = π/2 + θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    v = c .+ r .* [cos.(θs) sin.(θs)]'  # SMatrix{2,K}: locations of vertices

    return Polygon(v)
end

Polygon{K}(c::AbstractVector{<:Number},  # [x, y]: center of regular polygon
           r::Number,  # radius: distance from center to vertices
           θ::Number=0.0  # angle of first vertex
           ) where {K} =
   Polygon{K}(SVector{2}(c), r, θ)


Base.:(==)(s1::Polygon, s2::Polygon) = s1.v==s2.v && s1.n==s2.n  # assume sorted v
Base.isapprox(s1::Polygon, s2::Polygon) = s1.v≈s2.v && s1.n≈s2.n  # assume sorted v
Base.hash(s::Polygon, h::UInt) = hash(s.v, hash(s.n, hash(:Polygon, h)))

function level(x::SVector{2,<:Number}, s::Polygon)
    c = mean(s.v, dims=Val(2))  # center of mass

    d = sum(s.n .* (x .- c), dims=Val(1))
    r = sum(s.n .* (s.v .- c), dims=Val(1))
    @assert all(r .> 0)

    return 1 - maximum(d ./ r)
end

function surfpt_nearby(x::SVector{2,<:Number}, s::Polygon{K}) where {K}
    # Calculate the signed distances from x to edge lines.
    ∆xe = sum(s.n .* (x .- s.v), dims=Val(1))[1,:]  # SVector{K}: values of equations of edge lines
    abs∆xe = abs.(∆xe)  # SVector{K}

    # Determine if x is outside of edges, inclusive.
    sz = abs.((-)(bounds(s)...))  # SVector{2}
    onbnd = abs∆xe .≤ rtol(maximum(sz))  # SVector{K}
    isout = (∆xe.>0) .| onbnd  # SVector{K}

    # For x inside the polygon, it is easy to find the closest surface point: we can simply
    # pick the closest edge and find its point closest to x.
    # For x outside the polygon, there are many cases depending on how many edges x lies
    # outside of.  However, x that is sufficiently close to the polygon should be outside of
    # at most two edges.  (Outside two edges near the vertices, and one edge elsewhere.)
    # Therefore, we will consider only cout ≤ 2 below.
    cout = count(isout)
    if cout == 2  # x is outside two edges
        # We could choose to find ind corresponding to the two nonnegative ∆xe, but such an
        # operation leads to allocations because it does not create an SVector.  Instead,
        # find the closest vertex directly, which is the surface point we are looking for
        # for cout = 2.
        ∆xv = x .- s.v
        l∆xv = hypot.(∆xv[1,:], ∆xv[2,:])
        imin = argmin(l∆xv)
        surf = s.v[:,imin]
        imin₋₁ = mod1(imin-1,K)

        if onbnd[imin] && onbnd[imin₋₁]  # x is very close to vertex imin
            nout = s.n[:,imin] + s.n[:,imin₋₁]
        else
            nout = x - s.v[:,imin]
        end
        nout = normalize(nout)
    else  # cout ≤ 1 or cout ≥ 3
        # For cout ≥ 3, x is too far from the polygon and the result could be inaccurate.
        # (This case is not warned about, because logging macros expand to try/catch blocks
        # that reverse-mode AD tools such as Mooncake cannot differentiate through.)
        # Choose the closest edge to x.
        # If cout = 0, all ∆xe are negative, so the largest ∆xe is the smallest in magnitude
        # and corresponds to the closest edge.
        # If cout = 1, all but one ∆xe are negative, so again the largest ∆xe is the only
        # nonnegative one and corresponds to the edge outside which x lies.
        # Even for cout = 3, this seems to be a reasonable choice from a simple geometric
        # argument.
        imax = argmax(∆xe)
        vmax, nmax = s.v[:,imax], s.n[:,imax]

        ∆x = (nmax⋅(vmax-x)) .* nmax
        surf = x + ∆x
        nout = nmax
    end

    return surf, nout
end

translate(s::Polygon{K,K2}, ∆::SVector{2,<:Number}) where {K,K2} = Polygon{K,K2}(s.v .+ ∆, s.n)

function bounds(s::Polygon)
    l = minimum(s.v, dims=Val(2))[:,1]
    u = maximum(s.v, dims=Val(2))[:,1]

    return (l, u)
end
