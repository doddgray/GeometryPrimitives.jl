module GeometryPrimitives

using StaticArrays
using LinearAlgebra
using Statistics: mean

export Shape, Shape1, Shape2, Shape3
export level, surfpt_nearby, normal, bounds, translate
export drawshape, drawshape!

abstract type Shape{N,N²} end # a solid geometric shape in N dimensions (N² = N*N is needed in some shapes, e.g., Cuboid)
const Shape1 = Shape{1,1}
const Shape2 = Shape{2,4}
const Shape3 = Shape{3,9}

Base.ndims(::Shape{N}) where {N} = N

# Plotting stubs: methods are added by the GeometryPrimitivesMakieExt package extension,
# which is loaded automatically when Makie is loaded alongside this package.
function drawshape end
function drawshape! end

# relative tolerance * x (assumed ≥ 0) for approximate comparisons,
# defined to square root of machine precision like Base.rtoldefault
rtol(x::Real) = (y = float(x); sqrt(eps(typeof(y))) * y)

# Floating-point element type used to store the parameters of a shape constructed from
# arguments with the given element types.  Keeping this generic (instead of forcing
# Float64) lets shapes carry dual/tracked number types used by automatic differentiation.
promote_eltype(Ts::Type...) = promote_type(map(float, Ts)...)

# The following functions return Any due to the limitations of Julia's dispatch system.
# Therefore, always call them with return type assertions.  See
# https://discourse.julialang.org/t/extending-base-in-type-stably/5341/12
# https://github.com/JuliaLang/julia/issues/23210
level(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = level(SVector{N}(x), s)
Base.in(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = level(x,s) ≥ 0
surfpt_nearby(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = surfpt_nearby(SVector{N}(x), s)
normal(x::AbstractVector{<:Real}, s::Shape) = surfpt_nearby(x, s)[2]  # outward direction even for x inside s
translate(s::Shape{N}, ∆::AbstractVector{<:Real}) where {N} = translate(s, SVector{N}(∆))
# Shapes are immutable (for the sake of automatic differentiation), so there is no generic
# mutating implementation of translate(s, ∆::SVector); each shape type defines its own
# method that constructs a translated copy.

function orthoaxes(n::SVector{3,<:Real})
    u_temp = abs(n[3]) < abs(n[1]) ? SVector(0,0,1) : SVector(1,0,0)
    v = normalize(n × u_temp)
    u = v × n

    return u, v
end

include("hyper/hyper.jl")
include("planar/planar.jl")
include("prism/prism.jl")
include("util/util.jl")

end # module
