export Polygon, PolygonalPrism
export regpoly, isosceles

#= Polygon (for a base shape) =#

# Assume the followings for the polygon represented by Polygon:
# - The polygon is convex.
# - The vertices are listed in the counter-clockwise order around the origin.

mutable struct Polygon{K,K2,D,T} <: Shape{2,4,D,T}  # K2 = 2K
    v::SMatrix{K,2,T,K2}  # vertices
    n::SMatrix{K,2,T,K2}  # direction normals to edges
	l::SVector{2,T}		  # lower x,y bounds
	u::SVector{2,T}		  # upper x,y bounds
	sz::SVector{2,T}	  # "size" in x and y: sz = l-u
	rbnd::T				  # minimum distance before point is considered to be on boundary, rmin = √eps(T) * max(sz.data...)
    data::D  # auxiliary data
    # Polygon{K,K2,D,T}(v,n,data) where {K,K2,D,T<:Real} = new(v,n,data)  # suppress default outer constructor
end
# mutable struct Polygon{K,K2,D} <: Shape{2,4,D}  # K2 = 2K
#     v::SMatrix{K,2,Float64}  # vertices
#     n::SMatrix{K,2,Float64}  # direction normals to edges
#     data::D  # auxiliary data
#     Polygon{K,K2,D}(v,n,data) where {K,K2,D} = new(v,n,data)  # suppress default outer constructor
# end

function l_bnds_poly(v::SMatrix{K,2,T}) where {K,T<:Real}
	@tullio (min) res[b] :=  v[a,b]
	return SVector(res)
end
function u_bnds_poly(v::SMatrix{K,2,T}) where {K,T<:Real}
	@tullio (max) res[b] :=  v[a,b]
	return SVector(res)
end
# function Polygon(v::AbstractMatrix{T}, data::D=nothing) where {D,T<:Real}
#     K = size(v,1);


function ∆ϕ(∆v::SMatrix{K,2,<:Real}) where K
	# Calculate the increases in angle between neighboring edges.
	∆z = ∆v[:,1] + im * ∆v[:,2]  # SVector{K}: edge directions as complex numbers
	icurr = ntuple(identity, Val(K-1))
	inext = ntuple(x->x+1, Val(K-1))
	return angle.(∆z[SVector(inext)] ./ ∆z[SVector(icurr)])  # angle returns value between -π and π
end

function sort_v_if_needed(v::SMatrix{K,2,T}) where {K,T}
	w = v .- mean(v, dims=1)  # v in center-of-mass coordinates
	ϕ = mod.(atan.(w[:,2], w[:,1]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
	if !issorted(ϕ)	# TODO: make sort_verts shuffling fn with AD rules, currently unsorted verts would break differentiability
		# Do this only when ϕ is not sorted, because the following uses allocations.
		ind = MVector{K}(sortperm(ϕ))  # sortperm(::SVector) currently returns Vector, not MVector
		v = v[ind,:]  # SVector{K}: sorted v
	end
end

function n_norm(dv)
	local one_mone = [1, -1]
	@tullio n[i,j] := dv[i,3-j] * $one_mone[j] / sqrt( dv[i,1]^2 + dv[i,2]^2 ) nograd=one_mone
end

n_norm_fwd(dv) = Zygote.forwarddiff(n_norm,dv)

ChainRulesCore.@non_differentiable ∆ϕ(::Any)

function Polygon(v::SMatrix{K,2,T}, data::D=nothing) where {K,D,T<:Real}
    # Sort the vertices in the counter-clockwise direction

    # w = v .- mean(v, dims=1)  # v in center-of-mass coordinates
    # ϕ = mod.(atan.(w[:,2], w[:,1]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
    # if !issorted(ϕ)	# TODO: make sort_verts shuffling fn with AD rules, currently unsorted verts would break differentiability
    #     # Do this only when ϕ is not sorted, because the following uses allocations.
    #     ind = MVector{K}(sortperm(ϕ))  # sortperm(::SVector) currently returns Vector, not MVector
    #     v = v[ind,:]  # SVector{K}: sorted v
    # end
    ∆v = v - circshift(v,1)
    # Check all the angle increases are positive.  If they aren't, the polygon is not convex.
    # all(∆ϕ(∆v) .> 0) || throw("v = $v should represent vertices of convex polygon.")

	# outward normal directions to edges
	# one_mone = [	0. 1.
	# 				-1. 0. ]# SMatrix{2,2}(0., 1., -1., 0.)
    # n0 = ∆v * one_mone #[∆v[:,2] -∆v[:,1]]  # outward normal directions to edges
	# @tullio nnorm[k] := ∆v[k,j]^2 |> sqrt
	# @tullio n[k,j] := n0[k,j] / nnorm[k] # n = n0 ./ hypot.(n0[:,1],n0[:,2])  # normalize
	# n0 = [∆v[:,2] -∆v[:,1]]  # outward normal directions to edges
    # n = n0 ./ hypot.(n0[:,1],n0[:,2])  # normalize

	# n = n_norm(∆v)
	n = n_norm_fwd(∆v)
	l = l_bnds_poly(v)
	u = u_bnds_poly(v)
	sz = u-l
	rbnd = Zygote.@ignore(max(sz.data[1],sz.data[2])*Base.rtoldefault(T))
	# TODO: Ask why `bounds(::Polygon)` was prev. recomputed each time and why
	# the "size" (`sz`) in find_surfpt_nearby(Polygon was computed as abs.(l-u)?
	# Hopefully this doesn't break something or hurt performance
    return Polygon{K,2K,D,T}(v,SMatrix{K,2,T}(n),l,u,sz,rbnd,data) # Polygon{K,2K,D}(v,n,data)
end

function _∆xe_poly(x::SVector{2},v::SMatrix{K,2},n::SMatrix{K,2})::SVector{K}  where {K} #,T<:Real}
	@tullio out[k] := n[k,a] * ( x[a] - v[k,a] ) # edge line eqns for a K-point Polygon{K} `s`
end
function ∆xe_poly(x::AbstractVector{<:Real},s::Polygon{K})::SVector{K} where K
	_∆xe_poly(x,s.v,s.n)
end
bounds(s::Polygon) = (s.l, s.u)

# function onbnd(Δxe::SVector{K},s::Polygon{K})::SVector{K,Bool} where K 		# fixed wrong Δ, keeping in case this breaks
function onbnd(∆xe::SVector{K},s::Polygon{K})::SVector{K,Bool} where K
	map(x->abs(x)≤s.rbnd,∆xe)
end
function onbnd(x::SVector{2},s::Polygon)::SVector{K,Bool} where K
	map(x->abs(x)≤s.rbnd,∆xe_poly(x,s))
end
function cout(x::SVector{2},s::Polygon)::Int
	mapreduce((a,b)->Int(a|b),+,onbnd(x,s),map(isposdef,∆xe_poly(x,s)))
end
function cout(∆xe::SVector{K,<:Real},obd::SVector{K,Bool})::Int where K
	mapreduce((a,b)->Int(a|b),+,obd,map(isposdef,∆xe))
end
function imin2(x::SVector{2},s::Polygon{K})::Tuple{Int,Int} where K
	∆xv = x' .- s.v
	l∆xv = hypot.(∆xv[:,1], ∆xv[:,2])
	imin = argmin(l∆xv)
	return  imin ,  mod1(imin-1,K) # imin,imin₋₁
end

function _∆x_poly(x::SVector{2},v::SMatrix{K,2},n::SMatrix{K,2},kmax::Int)::SVector{2}  where {K} #,T<:Real}
	# @tullio Δx[i] := n[$kmax,a] * ( v[$kmax,a] - x[a] * n[$kmax,i]	# works but gradient doesn't vectorize
	( @tullio Δx_factor := n[$kmax,a] * ( v[$kmax,a] - x[a] ) )  * SVector(n[kmax,1],n[kmax,2])
end
∆x_poly(x::AbstractVector{<:Real},s::Polygon)::SVector{2} = _∆x_poly(x,s.v,s.n,argmax(_∆xe_poly(x,s.v,s.n)))

ChainRulesCore.@non_differentiable cout(::Any,::Any)
ChainRulesCore.@non_differentiable onbnd(::Any,::Any)
ChainRulesCore.@non_differentiable imin2(::Any,::Any)


# Polygon(v::AbstractMatrix{<:Real}, data=nothing) = (K = size(v,1); Polygon(SMatrix{K,2}(v), data))

Base.:(==)(s1::Polygon, s2::Polygon) = s1.v==s2.v && s1.n==s2.n && s1.data==s2.data  # assume sorted v
Base.isapprox(s1::Polygon, s2::Polygon) = s1.v≈s2.v && s1.n≈s2.n && s1.data==s2.data  # assume sorted v
Base.hash(s::Polygon, h::UInt) = hash(s.v, hash(s.n, hash(s.data, hash(:Polygon, h))))

Base.in(x::SVector{2,<:Real}, s::Polygon) = all(sum(s.n .* (x' .- s.v), dims=Val(2)) .≤ 0)
# Base.in(x::SVector{3,<:Real}, s::Polygon) = ( x2 = SVector(x[1], x[2]); all(sum(s.n .* (x2' .- s.v), dims=Val(2)) .≤ 0) )

function surfpt_nearby(x::SVector{2,T}, s::Polygon{K}) where {K,T<:Real}
    # Calculate the signed distances from x to edge lines.
	# ∆xe = ∆xe_poly(x,s)
	∆xe = _∆xe_poly(x,s.v,s.n)
	# Determine if x is outside of edges, inclusive.
	obd = Zygote.@ignore(onbnd(∆xe,s))
    # For x inside the polygon, it is easy to find the closest surface point: we can simply
    # pick the closest edge and find its point closest to x.
    # For x outside the polygon, there are many cases depending on how many edges x lies
    # outside of.  However, x that is sufficiently close to the polygon should be outside of
    # at most two edges.  (Outside two edges near the vertices, and one edge elsewhere.)
    # Therefore, we will consider only cout ≤ 2 below.
    co = Zygote.@ignore(cout(∆xe,obd))
    if co == 2  # x is outside two edges
        # We could choose to find ind corresponding to the two nonnegative ∆xe, but such an
        # operation leads to allocations because it does not create an SVector.  Instead,
        # find the closest vertex directly, which is the surface point we are looking for
        # for cout = 2.
        imin, imin₋₁ = imin2(x,s)
		surf = SVector{2}(s.v[imin,1],s.v[imin,2])
		@inbounds if obd[imin] && obd[imin₋₁]  # x is very close to vertex imin
            # nout = reinterpret(SVector{2,T}, normalize( [ s.n[imin,1]+s.n[imin₋₁,1],s.v[imin,2]+s.n[imin₋₁,2] ]))  #  s.n[imin,:] + s.n[imin₋₁,:]
			@inbounds nout = normalize( [ s.n[imin,1]+s.n[imin₋₁,1],s.n[imin,2]+s.n[imin₋₁,2] ])
		else
            # nout = reinterpret(SVector{2,T}, normalize( [ x[1] - s.v[imin,1] , x[2] - s.v[imin,2] ] ) )[1]
			@inbounds nout = normalize( [ x[1] - s.v[imin,1] , x[2] - s.v[imin,2] ] )
        end
    else  # cout ≤ 1 or cout ≥ 3
        # co ≤ 1 || @warn "x = $x is outside $cout edges: too far from polygon with vertices $(s.v); " *
        #                     "result could be inaccurate."
        # Choose the closest edge to x.
        # If cout = 0, all ∆xe are negative, so the largest ∆xe is the smallest in magnitude
        # and corresponds to the closest edge.
        # If cout = 1, all but one ∆xe are negative, so again the largest ∆xe is the only
        # nonnegative one and corresponds to the edge outside which x lies.
        # Even for cout = 3, this seems to be a reasonable choice from a simple geometric
        # argument.
        imax::Int = Zygote.@ignore(argmax(∆xe))
        surf = x + _∆x_poly(x,s.v,s.n,imax) # ∆x = (nmax⋅(vmax-x)) .* nmax
        nout = @inbounds SVector(s.n[imax,1],s.n[imax,2] )
    end

    return surf, nout
end

translate(s::Polygon{K,K2,D}, ∆::SVector{2,T}) where {K,K2,D,T<:Real} = Polygon{K,K2,D,T}(s.v .+ transpose(∆), s.n, s.data)





#= Factory methods =#
# Regular polygon
function regpoly(::Val{K},  # number of vertices
                 r::Real,  # distance between center and each vertex
                 θ::Real=π/2,  # angle from +x-direction towards first vertex; π/2 corresponds to +y-direction
                 c::SVector{2,T}=SVector(0.0,0.0),  # center location
                 data=nothing) where {K,T<:Real}
    ∆θ = 2π / K

    θs = θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    v = c' .+ r .* [cos.(θs) sin.(θs)]  # SMatrix{K,2}: locations of vertices

    return Polygon(v, data)
end

regpoly(k::Integer,  # number of vertices
        r::Real,  # radius: distance from center to vertices
        θ::Real=π/2,  # angle of first vertex
        c::AbstractVector{<:Real}=[0.0,0.0],  # [x, y]: center of regular polygon
        data=nothing) =
   regpoly(Val(k), r, θ, SVector{2}(c), data)

# Isosceles triangle
function isosceles(base::NTuple{2,SVector{2,<:Real}},
                   h::Real,
                   data=nothing)
    m = (base[1] + base[2]) / 2  # midpoint of base
    bvec = normalize(base[2] - base[1])  # unit direction of base
    hvec = [-bvec[2], bvec[1]]  # unit direction of height
    p = m + h.*hvec  # apex

    v = [base[1] base[2] p]'  # vertices

    return Polygon(v, data)
end

isosceles(base::NTuple{2,AbstractVector{<:Real}},  # (end point 1, end point 2): two end points of base
          h::Real,  # height drawn normal to base; direction is such that base pt 1, base pt 2, apex are put in counter-clockwise order
          data=nothing) = isosceles(SVector{2}.(base), h, data)

# To-dos: parallegram, rhombus, ...


#= Polygonal prism =#
const PolygonalPrism{K,K2,T} = Prism{Polygon{K,K2,Nothing,T}}

# Below, if we called PolygonalPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Polygon{K,K2,Nothing}}(c, ...) because PolygonalPrism = Prism{Polygon{K,K2,Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of PolygonalPrism(c, ...).
PolygonalPrism(c::SVector{3,<:Real},
               v::SMatrix{K,2,<:Real},  # 2D coordinates of base vertices in projected prism coordinates
               h::Real=Inf,
               a::SVector{3,<:Real}=SVector(0.0,0.0,1.0),
               data=nothing) where {K} =
    (â = normalize(a); Prism(c, Polygon(v), h, [orthoaxes(â)... â], data))

PolygonalPrism(c::AbstractVector{<:Real},  # center of prism
               v::AbstractMatrix{<:Real},  # vertices of base polygon
               h::Real=Inf,  # height of prism
               a::AbstractVector{<:Real}=[0.0,0.0,1.0],  # axis direction of prism
               data=nothing) =
    (K = size(v,1); PolygonalPrism(SVector{3}(c), SMatrix{K,2}(v), h, SVector{3}(a), data))

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::PolygonalPrism{K}) where {K}
    ax = inv(s.p)  # prism axes: columns are not only unit vectors, but also orthogonal
    v = [s.b.v @SVector(zeros(K))]  # SMatrix{K,3}: 3D vectices in prism axis coordinates
    w = v * ax'  # SMatrix{K,3}: vertices in external coordinates (equivalent to (ax * v')')

    return minimum(w, dims=Val(1))[1,:], maximum(w, dims=Val(1))[1,:]
end
