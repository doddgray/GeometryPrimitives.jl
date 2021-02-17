export Box

# Below, the box axes describe the directions of edges.  For example, if a₁, a₂, a₃ are the
# axes of a 3D box, then the volume inside the box is the collection of points x₁a₁ + x₂a₂
# + x₃a₃, where -rₙ ≤ xₙ ≤ rₙ.  This means aₙ and aₘ with n ≠ m span a face of the box.
#
# The box projection matrix (Box.p) is the matrix that produces the coordinates (x₁, x₂, x₃)
# of a given point.  In other words, for a point x inside the box, the nth entry of Box.p * x
# must have magnitude ≤ rₙ.  This means that Box.p * aₙ = eₙ, because aₙ has only the aₙ-component
# and the component must be 1.  This also means that each row of Box.p is orthogononal to
# two box axes, and therefore normal to the face spanned by the two box axes.  (Note that
# the rows of Box.p are not unit normals to the faces, becuase they are not unit vectors.)
mutable struct Box{N,N²,D,T} <: Shape{N,N²,D,T}
    c::SVector{N,T}  # center of box
    r::SVector{N,T}  # "radii" (semi-axes) in axis directions
    p::SMatrix{N,N,T,N²}  # projection matrix to box coordinates
	l::SVector{N,T}	# lower bounds
	u::SVector{N,T}	# upper bounds
    data::D  # auxiliary data
    #Box{N,N²,D,T}(c,r,p,data) where {N,N²,D,T<:Real} = new(c,r,p,data)  # suppress default outer constructor
end

signmatrix(b::Box{1}) = SMatrix{1,1}(1)
signmatrix(b::Box{2}) = SMatrix{2,2}(1,1, -1,1)
signmatrix(b::Box{3}) = SMatrix{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)
signmatrix(v::SVector{1}) = SMatrix{1,1}(1)
signmatrix(v::SVector{2}) = SMatrix{2,2}(1,1, -1,1)
signmatrix(v::SVector{3}) = SMatrix{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)
@non_differentiable signmatrix(::Any)

function Box(c::SVector{N,T},
		    d::SVector{N,T},
		    axes::SMatrix{N,N,<:Real}=SMatrix{N,N,T}(I),
		    data::D=nothing) where {N,D,T<:Real}
	r = 0.5d
	@tullio axnorm[i] := axes[i,j]^2 |> sqrt
	@tullio p_inv[i,j] := axes[i,j] / axnorm[j]
	# p_inv = axes ./ sqrt.(sum(abs2,axes,dims=1))
	smr = Zygote.@ignore(signmatrix(r))
	@tullio (max) m[j] := p_inv[i,k] * r[i] * smr[i,j] nograd=smr
    Box{N,N*N,D,T}(c, 0.5d, SMatrix(inv(p_inv)),c-SVector(m),c+SVector(m),data)
end

Box(c::AbstractVector{T},  # center of box
    d::AbstractVector{T},  # size of box in axis directions
    axes::AbstractMatrix{<:Real}=Matrix{T}(I,length(c),length(c)),  # columns are axes vectors (each being parallel to two sets of faces in 3D)
    data=nothing) where T<:Real =
    (N = length(c); Box(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Box, b2::Box) = b1.c==b2.c && b1.r==b2.r && b1.p==b2.p && b1.data==b2.data
Base.isapprox(b1::Box, b2::Box) = b1.c≈b2.c && b1.r≈b2.r && b1.p≈b2.p && b1.data==b2.data
Base.hash(b::Box, h::UInt) = hash(b.c, hash(b.r, hash(b.p, hash(b.data, hash(:Box, h)))))

function bounds(b::Box)
    # A = inv(b.p) .* b.r'
    # m = maximum(abs.(A * signmatrix(b)), dims=2)[:,1] # extrema of all 2^N corners of the box
    # # m = maximum(abs.(A * signmatrix(b)), dims=2)[:,1] # extrema of all 2^N corners of the box
    # # m = maximum(abs.(Array((inv(b.p) .* b.r') * signmatrix(b))), dims=2)[:,1] # extrema of all 2^N corners of the box
	#
	# # A = ( inv(b.p) .* b.r' ) * signmatrix(b)
    # # m = maximum(abs.(A * signmatrix(b)), dims=2)[:,1] # extrema of all 2^N corners of the box
    # return (b.c-m,b.c+m)
	return b.l, b.u
end

function Base.in(x::SVector{N,<:Real}, b::Box{N}) where {N}
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false  # boundary is considered inside
    end
    return true
end

function f_onbnd(bin::Box{N,N²,D,T},absdin) where {N,N²,D,T<:Real}
	b = Zygote.@ignore(bin)
	absd = Zygote.@ignore(absdin)
	abs.(b.r.-absd) .≤ Base.rtoldefault(T) .* b.r  # basically b.r .≈ absd but faster
end

f_isout(b,absd) =  (isout = Zygote.@ignore ((b.r.<absd) .| f_onbnd(b,absd) ); isout)
f_signs(d) =  (signs = Zygote.@ignore(sign.(d)'); signs) #(copysign.(1.0,d)'); signs)

@non_differentiable f_onbnd(::Any)
@non_differentiable f_isout(::Any)
@non_differentiable f_signs(::Any)

function surfpt_nearby(x::SVector{N,T}, b::Box{N}) where {N,T<:Real}
    ax = inv(b.p)  # axes: columns are unit vectors
	n0 = @inbounds b.p ./  [ sqrt(b.p[1,1]^2 + b.p[1,2]^2) sqrt(b.p[2,1]^2 + b.p[2,2]^2)  ]
    d= Array(b.p * (x - b.c))
    cosθ = diag(n0*ax)
	n = n0 .* f_signs(d)
    # n = n0 .* copysign.(1.0,d)' #f_signs(d)
    absd = abs.(d)
    ∆ = (b.r .- absd) .* cosθ
    onbnd = f_onbnd(b,absd)
    isout = f_isout(b,absd)
    # projbnd =  all(.!isout .| onbnd)
    if count(isout) == 0  # x strictly inside box; ∆ all positive
        l∆x, i = findmin(∆)  # find closest face
        nout = n[i,:]
        ∆x = l∆x * nout
    else  # x outside box or on boundary in one or multiple directions
        ∆x = n' * (∆ .* isout)  # project out .!isout directions

        # Below, (.!isout .| onbnd) tests if each entry of ∆ (not ∆x) is either projected
        # out or on the boundary.  If the ith entry of ∆ is projected out, ∆x (not ∆) does
        # not have any component in the n[i,:] direction.  Even if the ith entry of ∆ is not
        # projected out, if it is too small then ∆x barely has any component in the n[i,:]
        # direction.  If all the entries of ∆ satisfy one of the two conditions, ∆x ≈ 0 and
        # we cannot use -∆x for nout.  In that case, take n'*onbnd as nout.  When x is
        # outside the box only in one dimension, n'*onbnd is basically the row of n along
        # that dimension, which reduces to normalize(-∆x) even if ∆x is very small.
        # nout = all(.!isout .| onbnd) ? n'*onbnd : -∆x
        # nout = normalize(nout)
		if all(.!isout .| onbnd)
            nout0 = n' * onbnd
        else
            nout0 = -∆x
        end
        nout = nout0 / norm(nout0)
    end

    return x+∆x, nout
end

translate(b::Box{N,N²,D}, ∆::SVector{N,T}) where {N,N²,D,T<:Real} = Box{N,N²,D,T}(b.c+∆, b.r, b.p, b.data)
