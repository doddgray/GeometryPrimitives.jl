# Gradient (automatic differentiation) tests for GeometryPrimitives.
#
# The differentiable API surface is:
#   - level(x, shape)            w.r.t. x and w.r.t. shape parameters
#   - surfpt_nearby(x, shape)    w.r.t. x and w.r.t. shape parameters
#   - volfrac(vxl, nout, r₀)     w.r.t. the voxel corners and the plane parameters
#
# Gradients are computed with Enzyme.jl (forward and reverse mode) and Mooncake.jl
# (reverse mode) through DifferentiationInterface.jl, and verified against central finite
# differences (FiniteDifferences.jl).

using GeometryPrimitives
using StaticArrays
using LinearAlgebra
using Test
using Random: MersenneTwister

using DifferentiationInterface
using ADTypes: AutoEnzyme, AutoMooncake, AutoFiniteDifferences
import Enzyme, Mooncake, FiniteDifferences

const GRAD_RNG = MersenneTwister(42)

const FD = AutoFiniteDifferences(; fdm=FiniteDifferences.central_fdm(5, 1))
const AD_BACKENDS = (
    "Enzyme reverse" => AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse),
                                     function_annotation=Enzyme.Const),
    "Enzyme forward" => AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward),
                                     function_annotation=Enzyme.Const),
    "Mooncake" => AutoMooncake(; config=nothing),
)

# Test the gradient of f at p against finite differences for every AD backend.
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

# Scalar objectives
flevel(sh)  = x -> level(SVector{length(x)}(x), sh)
fsurfpt(sh) = x -> sum(sum, surfpt_nearby(SVector{length(x)}(x), sh))
function fvolfrac(sh, ::Val{N}) where {N}
    return function (x)
        xs = SVector{N}(x)
        δ = ones(SVector{N,eltype(x)})
        r₀, nout = surfpt_nearby(xs, sh)
        return volfrac((xs - δ, xs + δ), nout, r₀)
    end
end

@testset verbose = true "AD gradients" begin
    @testset "x-gradients, 2D shapes" begin
        @testset "$name" for (name, make, np) in SHAPES2D
            sh = make(rand(GRAD_RNG, np))
            x = rand(GRAD_RNG, 2)
            test_grad(flevel(sh), x)
            test_grad(fsurfpt(sh), x)
            test_grad(fvolfrac(sh, Val(2)), x)
        end
    end

    @testset "x-gradients, 3D shapes" begin
        @testset "$name" for (name, make, np) in SHAPES3D
            sh = make(rand(GRAD_RNG, np))
            x = rand(GRAD_RNG, 3)
            test_grad(flevel(sh), x)
            test_grad(fsurfpt(sh), x)
            test_grad(fvolfrac(sh, Val(3)), x)
        end
    end

    @testset "shape-parameter gradients, 2D shapes" begin
        @testset "$name" for (name, make, np) in SHAPES2D
            x₀ = SVector{2}(rand(GRAD_RNG, 2))
            p = rand(GRAD_RNG, np)
            test_grad(p -> level(x₀, make(p)), p)
            test_grad(p -> sum(sum, surfpt_nearby(x₀, make(p))), p)
            test_grad(p -> begin
                r₀, nout = surfpt_nearby(x₀, make(p))
                δ = ones(SVector{2,eltype(p)})
                volfrac((x₀ - δ, x₀ + δ), nout, r₀)
            end, p)
        end
    end

    @testset "shape-parameter gradients, 3D shapes" begin
        @testset "$name" for (name, make, np) in SHAPES3D
            x₀ = SVector{3}(rand(GRAD_RNG, 3))
            p = rand(GRAD_RNG, np)
            test_grad(p -> level(x₀, make(p)), p)
            test_grad(p -> sum(sum, surfpt_nearby(x₀, make(p))), p)
            test_grad(p -> begin
                r₀, nout = surfpt_nearby(x₀, make(p))
                δ = ones(SVector{3,eltype(p)})
                volfrac((x₀ - δ, x₀ + δ), nout, r₀)
            end, p)
        end
    end
end
