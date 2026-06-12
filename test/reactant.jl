# Optional Reactant.jl tests: compile level-set evaluations and their Enzyme gradients to
# XLA via Reactant and compare against finite differences.
#
# These tests are not part of the default test suite because Reactant is a heavy
# dependency (it downloads an XLA runtime) and only supports a subset of platforms.
# Run them with:
#
#   GP_TEST_REACTANT=true julia --project=<env with Reactant, Enzyme, FiniteDifferences> test/reactant.jl
#
# or set GP_TEST_REACTANT=true when running the package test suite in such an environment.
#
# Notes on coverage: Reactant traces Julia code with TracedRNumber values, so only
# branch-free code paths can be compiled.  level() for Ball, Cuboid, Ellipsoid, and Polygon
# is branch-free (max/abs/sqrt reductions) and works; surfpt_nearby() and volfrac() branch
# on values and are not Reactant-traceable.  The type-parametric shape structs are what
# make this possible: shapes constructed inside a traced function carry TracedRNumber
# parameters.

using GeometryPrimitives
using StaticArrays
using Test
import Reactant, Enzyme, FiniteDifferences

# The objectives below index the traced parameter vector (p[1], p[2], ...) to build
# StaticArrays; each such getindex creates a TracedRNumber, which is exactly what we want
# during tracing, so allow it.
Reactant.allowscalar(true)

# Scalar objectives p::AbstractVector -> Real, differentiable w.r.t. p.
# Each constructs a shape from p (exercising traced shape construction) or evaluates at a
# traced point.
const x₀3 = SVector(0.3, 0.4, 0.5)
const x₀2 = SVector(0.3, 0.4)
rotmat2(θ) = SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ))

f_ball_x(x)       = level(SVector(x[1],x[2],x[3]), Ball(SVector(0.1,0.2,0.3), 1.1))
f_cuboid_x(x)     = level(SVector(x[1],x[2],x[3]), Cuboid(SVector(0.1,0.2,0.3), SVector(1.0,1.2,1.4)))
f_cyl_x(x)        = level(SVector(x[1],x[2],x[3]), Cylinder(SVector(0.1,0.2,0.3), 0.8, 1.7))
f_ball_p(p)       = level(x₀3, Ball(SVector(p[1],p[2],p[3]), p[4]))
f_cuboid_p(p)     = level(x₀3, Cuboid(SVector(p[1],p[2],p[3]), SVector(p[4],p[5],p[6])))
f_ellipsoid_p(p)  = level(x₀2, Ellipsoid(SVector(p[1],p[2]), SVector(p[3],p[4]), rotmat2(p[5])))

const REACTANT_CASES = (
    ("level Ball3 / x",           f_ball_x,      [0.5, 0.6, 0.7]),
    ("level Cuboid3 / x",         f_cuboid_x,    [0.5, 0.6, 0.7]),
    ("level Cylinder / x",        f_cyl_x,       [0.5, 0.6, 0.7]),
    ("level Ball3 / params",      f_ball_p,      [0.1, 0.2, 0.3, 1.1]),
    ("level Cuboid3 / params",    f_cuboid_p,    [0.1, 0.2, 0.3, 1.0, 1.2, 1.4]),
    ("level Ellipsoid2 / params", f_ellipsoid_p, [0.1, 0.2, 1.0, 1.4, 0.3]),
)

@testset verbose = true "Reactant + Enzyme gradients" begin
    fdm = FiniteDifferences.central_fdm(5, 1)
    @testset "$name" for (name, f, p) in REACTANT_CASES
        p_ra = Reactant.to_rarray(p)

        # Compiled primal must match the plain Julia evaluation.
        f_compiled = Reactant.@compile f(p_ra)
        @test Float64(f_compiled(p_ra)) ≈ f(p) rtol = 1e-12

        # Compiled Enzyme reverse-mode gradient must match finite differences.
        gradf(q) = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), q)[1]
        grad_compiled = Reactant.@compile gradf(p_ra)
        g = Array(grad_compiled(p_ra))
        g_fd = FiniteDifferences.grad(fdm, f, p)[1]
        @test g ≈ g_fd rtol = 1e-6
    end
end
