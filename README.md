# GeometryPrimitives

[![CI](https://github.com/stevengj/GeometryPrimitives.jl/workflows/CI/badge.svg)](https://github.com/stevengj/GeometryPrimitives.jl/actions)
[![Codecov](http://codecov.io/github/stevengj/GeometryPrimitives.jl/coverage.svg?branch=master)](http://codecov.io/github/stevengj/GeometryPrimitives.jl?branch=master)

This package provides a set of geometric primitive types (balls, cuboids, cylinders, and
so on) and operations on them designed to enable piecewise definition of functions,
especially for finite-difference and finite-element simulations, in the Julia language.

For example, suppose that you are discretizing a PDE like the Poisson equation ∇⋅c∇u = f,
and you want to provide a simple user interface for the user to specify the function `c(x)`.
In many applications, `c` will be piecewise constant, and you want to be able to specify
`c = 1` in one box, `c = 2` in some cylinders, etcetera.   The GeometryPrimitives package
allows the user to provide a list of shapes with associated data (in this case, the value of
`c`) to define such a `c(x)`.

Furthermore, the application to discretized simulations imposes a couple of additional
requirements:

* One needs to be able to evaluate `c(x)` a huge number of times (once for every point on a
grid).  So, we provide a fast O(log n) K-D tree data structure for rapid searching of shapes.

* Often, one wants to compute the *average* of `c(x)` over a voxel, so we provide routines
for rapid *approximate* voxel averages.

* Often, one needs not only the value `c(x)` but the normal vector to the nearest shape, so
we provide normal-vector computation.

This package was inspired by the geometry utilities in Steven G. Johnson's [Libctl]
(http://ab-initio.mit.edu/wiki/index.php/Libctl) package.

## Automatic differentiation

The shape types are immutable and parametrized on their numeric element type, and the
differentiable parts of the API contain no hidden `Float64` conversions, so the package
composes with Julia automatic-differentiation tools.  The differentiable surface is:

* `level(x, shape)` — w.r.t. the query point `x` and the shape parameters,
* `surfpt_nearby(x, shape)` (and `normal`) — w.r.t. `x` and the shape parameters,
* `volfrac(vxl, nout, r₀)` — w.r.t. the voxel corners and the cutting-plane parameters,

where "w.r.t. shape parameters" means differentiating through the shape constructors
(`Ball`, `Cuboid`, `Ellipsoid`, `Polygon`, `Sector`, `Cylinder`, `PolygonalPrism`,
`SectoralPrism`, …).

Gradients are tested against finite differences with:

* [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) (forward and reverse mode),
* [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) (reverse mode),

both used directly or through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
(see `test/grads.jl`), and with
[Reactant.jl](https://github.com/EnzymeAD/Reactant.jl) for XLA-compiled gradients of the
branch-free `level` functions (see `test/reactant.jl`; run with `GP_TEST_REACTANT=true`).

These AD backends have a high first-call compilation latency for this
`StaticArrays`-heavy code, so the gradient tests are split into four groups —
`x2d`, `x3d`, `param2d`, `param3d` (query-point vs. shape-parameter gradients, in 2D
and 3D) — selectable through the `GP_GRAD_GROUPS` environment variable so that each can be
run in its own process to keep memory bounded:

```sh
GP_GRAD_GROUPS=x2d     julia --project=test test/grads.jl
GP_GRAD_GROUPS=param3d julia --project=test test/grads.jl
```

The default package test suite (`julia --project -e 'using Pkg; Pkg.test()'`) runs the two
2D groups; set `GP_GRAD_GROUPS=all` (or `GP_TEST_AD_FULL=true`) to also exercise the 3D
groups.  The backend set can be narrowed similarly with `GP_GRAD_BACKENDS` (a subset of
`enzyme_reverse,enzyme_forward,mooncake`).  Benchmarks of primal vs. gradient evaluation
across the backends live in `benchmark/benchmarks.jl`.

Note that `surfpt_nearby` and `volfrac` select branches (nearest face, voxel/plane
crossing cases, …) based on the input values; their outputs are continuous and piecewise
differentiable, and AD returns the gradient of the active branch.

## Plotting

Loading [Makie](https://docs.makie.org) together with this package activates a package
extension that provides `drawshape(shape)` / `drawshape!(shape)` recipes for 2D shapes
(and 2D cross-sections of 3D shapes via `CrossSection`).
