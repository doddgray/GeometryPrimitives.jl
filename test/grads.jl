# Gradient (automatic differentiation) tests for GeometryPrimitives.
#
# The differentiable API surface is:
#   - level(x, shape)            w.r.t. x and w.r.t. shape parameters
#   - surfpt_nearby(x, shape)    w.r.t. x and w.r.t. shape parameters
#   - volfrac(vxl, nout, r₀)     w.r.t. the query point (via surfpt_nearby) and shape params
#
# Gradients are computed with Enzyme.jl (forward and reverse mode) and Mooncake.jl
# (reverse mode) through DifferentiationInterface.jl, and verified against central finite
# differences (FiniteDifferences.jl).
#
# These AD backends have a high first-call compilation latency for this StaticArrays-heavy
# code (tens of seconds per (shape, function, backend) combination), and the compiled-code
# caches accumulate within a process.  The tests are therefore organized into four groups
#
#     x2d  x3d  param2d  param3d
#
# selectable through the GP_GRAD_GROUPS environment variable (comma-separated; default
# "all").  Running each group in its own Julia process keeps memory bounded and lets CI
# parallelize, e.g.
#
#     GP_GRAD_GROUPS=x2d     julia --project test/grads.jl
#     GP_GRAD_GROUPS=param3d julia --project test/grads.jl
#
# The backend set can likewise be narrowed with GP_GRAD_BACKENDS (comma-separated subset of
# "enzyme_reverse", "enzyme_forward", "mooncake"; default all three).

using GeometryPrimitives
using StaticArrays
using LinearAlgebra
using Test
using Random: MersenneTwister

using DifferentiationInterface
using ADTypes: AutoEnzyme, AutoMooncake, AutoFiniteDifferences
import Enzyme, Mooncake, FiniteDifferences

const FD = AutoFiniteDifferences(; fdm=FiniteDifferences.central_fdm(5, 1))

const ALL_BACKENDS = (
    "enzyme_reverse" => ("Enzyme reverse",
                         AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse),
                                      function_annotation=Enzyme.Const)),
    "enzyme_forward" => ("Enzyme forward",
                         AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward),
                                      function_annotation=Enzyme.Const)),
    "mooncake" => ("Mooncake", AutoMooncake(; config=nothing)),
)

# Resolve which backends to run from GP_GRAD_BACKENDS (default: all).
function selected_backends()
    want = get(ENV, "GP_GRAD_BACKENDS", "all")
    keys = want == "all" ? first.(ALL_BACKENDS) : split(want, ',')
    return [name_be for (key, name_be) in ALL_BACKENDS if key in keys]
end
const AD_BACKENDS = selected_backends()

# Test the gradient of f at p against finite differences for every selected AD backend.
function test_grad(f, p; rtol=1e-5, atol=1e-8)
    g_fd = DifferentiationInterface.gradient(f, FD, p)
    @testset "$name" for (name, backend) in AD_BACKENDS
        g = DifferentiationInterface.gradient(f, backend, p)
        @test isapprox(g, g_fd; rtol=rtol, atol=atol)
    end
    return nothing
end

# Vertices of a regular K-gon inscribed in the unit circle (offset angle keeps vertices
# away from the ±x/±y axes, where vertex sorting could change order under perturbation).
function circpts(K)
    ϕ = range(0.05, 1.95, length=K+1)
    return SMatrix{2,K}(ntuple(i -> isodd(i) ? cospi(ϕ[(i+1)÷2]) : sinpi(ϕ[i÷2]), Val(2K))...)
end
rotmat2(θ) = SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ))
normcols(A::SMatrix{N,N}) where {N} = A ./ sqrt.(sum(abs2, A, dims=Val(1)))

#= Shape factory functions: param vector p -> shape =#
# 2D shapes
make_ball2(p)      = Ball(SVector(p[1],p[2]), p[3])
make_cuboid2(p)    = Cuboid(SVector(p[1],p[2]), SVector(p[3],p[4]), normcols(SMatrix{2,2}(p[5],p[6],p[7],p[8])))
make_ellipsoid2(p) = Ellipsoid(SVector(p[1],p[2]), SVector(p[3],p[4]), rotmat2(p[5]))
make_pgon5(p)      = Polygon(0.1*SMatrix{2,5}(ntuple(k->p[k], Val(10))) + circpts(5))
make_pgon10(p)     = Polygon(0.05*SMatrix{2,10}(ntuple(k->p[k], Val(20))) + circpts(10))
make_regpoly8(p)   = Polygon{8}(SVector(p[1],p[2]), p[3], p[4])
make_sector(p)     = Sector(SVector(p[1],p[2]), p[3], p[4], p[5])

# 3D shapes
make_ball3(p)      = Ball(SVector(p[1],p[2],p[3]), p[4])
make_cuboid3(p)    = Cuboid(SVector(p[1],p[2],p[3]), SVector(p[4],p[5],p[6]),
                            normcols(SMatrix{3,3}(ntuple(k->p[6+k], Val(9)))))
make_ellipsoid3(p) = Ellipsoid(SVector(p[1],p[2],p[3]), SVector(p[4],p[5],p[6]))
make_cylinder(p)   = Cylinder(SVector(p[1],p[2],p[3]), p[4], p[5], SVector(p[6],p[7],p[8]))
make_polyprism(p)  = PolygonalPrism(SVector(p[13],p[14],p[15]),
                                    0.1*SMatrix{2,6}(ntuple(k->p[k], Val(12))) + circpts(6),
                                    p[16], SVector(p[17],p[18],p[19]))
make_sectprism(p)  = SectoralPrism(SVector(p[1],p[2],p[3]), p[4], p[5], p[6], p[7],
                                   SVector(p[8],p[9],p[10]))

const SHAPES2D = (
    ("Ball2",       make_ball2,      3),
    ("Cuboid2",     make_cuboid2,    8),
    ("Ellipsoid2",  make_ellipsoid2, 5),
    ("Polygon5",    make_pgon5,     10),
    ("Polygon10",   make_pgon10,    20),
    ("RegPoly8",    make_regpoly8,   4),
    ("Sector",      make_sector,     5),
)
const SHAPES3D = (
    ("Ball3",          make_ball3,       4),
    ("Cuboid3",        make_cuboid3,    15),
    ("Ellipsoid3",     make_ellipsoid3,  6),
    ("Cylinder",       make_cylinder,    8),
    ("PolygonalPrism", make_polyprism,  19),
    ("SectoralPrism",  make_sectprism,  10),
)

# Differentiating *through* the Polygon constructor (sortperm + convexity check) is
# correct for every backend but very expensive to compile with Enzyme reverse mode, and the
# cost grows with the vertex count.  The param-gradient groups therefore use a single small
# polygon (Polygon5) as the representative polygon-construction case and drop the larger
# Polygon10 and RegPoly8 (which are still covered by the x-gradient groups, where the
# constructor is not differentiated).
const PARAM_SHAPES2D = filter(s -> first(s) ∉ ("Polygon10", "RegPoly8"), SHAPES2D)

# δ for the voxel half-width in the volfrac objective: large enough that the voxel straddles
# the surface for the random query points used here.
δvec(::Val{N}, T) where {N} = ones(SVector{N,T})

#= Test groups =#
function run_x_gradients(shapes, ::Val{N}) where {N}
    rng = MersenneTwister(42)
    @testset "$name" for (name, make, np) in shapes
        sh = make(rand(rng, np))
        x = rand(rng, N)
        @testset "level"   test_grad(x -> level(SVector{N}(x), sh), x)
        @testset "surfpt"  test_grad(x -> sum(sum, surfpt_nearby(SVector{N}(x), sh)), x)
        @testset "volfrac" test_grad(x -> begin
            xs = SVector{N}(x)
            r₀, nout = surfpt_nearby(xs, sh)
            volfrac((xs - δvec(Val(N),eltype(x)), xs + δvec(Val(N),eltype(x))), nout, r₀)
        end, x)
    end
end

function run_param_gradients(shapes, ::Val{N}) where {N}
    rng = MersenneTwister(7)
    @testset "$name" for (name, make, np) in shapes
        x₀ = SVector{N}(rand(rng, N))
        p = rand(rng, np)
        @testset "level"   test_grad(p -> level(x₀, make(p)), p)
        @testset "surfpt"  test_grad(p -> sum(sum, surfpt_nearby(x₀, make(p))), p)
        @testset "volfrac" test_grad(p -> begin
            r₀, nout = surfpt_nearby(x₀, make(p))
            δ = δvec(Val(N), eltype(p))
            volfrac((x₀ - δ, x₀ + δ), nout, r₀)
        end, p)
    end
end

const GRAD_GROUPS = Dict(
    "x2d"     => () -> run_x_gradients(SHAPES2D, Val(2)),
    "x3d"     => () -> run_x_gradients(SHAPES3D, Val(3)),
    "param2d" => () -> run_param_gradients(PARAM_SHAPES2D, Val(2)),
    "param3d" => () -> run_param_gradients(SHAPES3D, Val(3)),
)
const GROUP_TITLES = Dict(
    "x2d"     => "x-gradients, 2D shapes",
    "x3d"     => "x-gradients, 3D shapes",
    "param2d" => "shape-parameter gradients, 2D shapes",
    "param3d" => "shape-parameter gradients, 3D shapes",
)
const GROUP_ORDER = ("x2d", "x3d", "param2d", "param3d")

function selected_groups()
    want = get(ENV, "GP_GRAD_GROUPS", "all")
    want == "all" && return GROUP_ORDER
    return Tuple(g for g in GROUP_ORDER if g in split(want, ','))
end

@testset verbose = true "AD gradients" begin
    @testset verbose = true "$(GROUP_TITLES[g])" for g in selected_groups()
        GRAD_GROUPS[g]()
    end
end
