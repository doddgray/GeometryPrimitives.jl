using GeometryPrimitives
using StaticArrays
using LinearAlgebra
using Random: MersenneTwister
using Statistics: mean
using Test

const one⁻ = 1 - 1e-8  # scale factor slightly less than 1
const one⁺ = 1 + 1e-8  # scare factor slightly greater than 1
const one⁻⁻, one⁺⁺ = 0.9, 1.1  # (scale factor less than 1, scale factor greater than 1)

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))
const rng = MersenneTwister(0) # test with reproducible pseudorandom numbers

randnb(lb::Real,ub::Real) = randn(rng)*(ub-lb) + 0.5*(lb+ub)
randnb(lb::AbstractVector,ub::AbstractVector) = map(randnb, lb, ub)
function inbounds(x,lb,ub)
    length(x) == length(lb) == length(lb) || throw(BoundsError())
    for i in eachindex(x)
        @inbounds lb[i] ≤ x[i] ≤ ub[i] || return false
    end
    return true
end

"check the bounding box of s with randomized trials"
function checkbounds(s::Shape{N}, ntrials=10^4) where {N}
    lb,ub = bounds(s)
    for i = 1:ntrials
        x = randnb(lb,ub)
        x ∉ s || inbounds(x,lb,ub) || return false  # return false if s - (bounding box) is nonempty
    end
    return true
end

function checktree(t::KDTree{N}, slist::Vector{<:Shape{N}}, ntrials=10^3) where {N}
    lb = SVector{N}(fill(Inf,N))
    ub = SVector{N}(fill(-Inf,N))
    for i in eachindex(slist)
        lbi,ubi = bounds(slist[i])
        lb = min.(lb,lbi)
        ub = max.(ub,ubi)
    end
    for i = 1:ntrials
        x = randnb(lb,ub)
        st = findfirst(x, t)
        sl = findfirst(x, slist)
        st == sl || return false
    end
    return true
end

@testset "GeometryPrimitives" begin

include("ball.jl")
include("cuboid.jl")
include("ellipsoid.jl")
include("cross_section.jl")
include("polygon.jl")
include("sector.jl")
include("cylinder.jl")
include("polygonal_prism.jl")
include("sectoral_prism.jl")
include("kdtree.jl")
include("periodize.jl")
include("vxlcut.jl")

# Automatic-differentiation gradient tests.  These have a high AD-compilation cost
# (differentiating through the shape constructors — e.g. the Polygon constructor — is
# especially heavy for the param2d/param3d groups), so:
#   - set GP_TEST_AD=false to skip them (CI runs all groups in a separate parallel `ad`
#     job, one process per group);
#   - by default only the x2d group runs as a quick smoke test; GP_TEST_AD_FULL=true (or
#     GP_GRAD_GROUPS=all) runs all four groups.  GP_GRAD_GROUPS, if already set, is
#     respected as-is.  See grads.jl.
if get(ENV, "GP_TEST_AD", "true") == "true"
    if !haskey(ENV, "GP_GRAD_GROUPS")
        ENV["GP_GRAD_GROUPS"] = get(ENV, "GP_TEST_AD_FULL", "false") == "true" ? "all" : "x2d"
    end
    include("grads.jl")
end

# Optional: Reactant tests (heavy XLA dependency; see test/reactant.jl for details).
if get(ENV, "GP_TEST_REACTANT", "false") == "true"
    include("reactant.jl")
end

# Optional: Makie visualization tests (heavy; require CairoMakie in the active environment,
# e.g. the `viz` project — see test/visualize.jl).  Run with GP_TEST_VIZ=true, or directly
# via `julia --project=viz test/visualize.jl`.
if get(ENV, "GP_TEST_VIZ", "false") == "true"
    include("visualize.jl")
end

end  # @testset "GeometryPrimitives"
