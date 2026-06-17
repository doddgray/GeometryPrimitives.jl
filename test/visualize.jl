# Optional visualization tests.  For every geometry primitive these render an example shape
# with the Makie extension and save it as a PNG, checking that the file is produced:
#
#   * 2D shapes (Ball, Cuboid, Ellipsoid, Polygon, Sector): a 2D plot.
#   * 3D shapes (Ball, Cuboid, Ellipsoid, Cylinder, PolygonalPrism, SectoralPrism): a 3D
#     perspective view plus axis-aligned 2D slice plots.
#   * combined scenes mixing several primitives, in 2D and 3D.
#
# These are heavy (they pull in Makie), so they are not part of the default test suite.
# Run them with CairoMakie available, e.g.
#
#     julia --project=viz test/visualize.jl
#
# or, inside the package test suite, with GP_TEST_VIZ=true (CairoMakie must be installed in
# the active environment).  PNGs are written to GP_VIZ_DIR (default test/viz_output).

using GeometryPrimitives
using StaticArrays
using LinearAlgebra
using Test
using CairoMakie
CairoMakie.activate!(type="png")

const VIZ_DIR = get(ENV, "GP_VIZ_DIR", joinpath(@__DIR__, "viz_output"))
mkpath(VIZ_DIR)

# Save `fig` to VIZ_DIR/name and assert the PNG was written; print the path so generated
# plots can be located (and surfaced) as they are produced.
function savefig(name, fig)
    path = joinpath(VIZ_DIR, name)
    save(path, fig)
    @test isfile(path) && filesize(path) > 0
    println("saved ", path)
    return path
end

rot2(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
function rotz(θ)
    c, s = cos(θ), sin(θ)
    return @SMatrix [c -s 0.0; s c 0.0; 0.0 0.0 1.0]
end

# Example shapes (kept convex so the level-set ray-casting in the extension is exact).
const SHAPES_2D = [
    ("ball",      Ball([0.0,0.0], 1.0),                                   :dodgerblue),
    ("cuboid",    Cuboid([0.2,0.0], [2.0,1.2], rot2(deg2rad(25))),        :tomato),
    ("ellipsoid", Ellipsoid([0.0,0.0], [2.0,1.0], rot2(deg2rad(30))),     :mediumpurple),
    ("polygon",   Polygon{6}([0.0,0.0], 1.2, deg2rad(10)),               :seagreen),
    ("sector",    Sector([0.0,0.0], 1.3, deg2rad(20), deg2rad(110)),      :darkorange),
]

const _pentagon = Polygon{5}([0.0,0.0], 0.9, 0.0).v  # SMatrix{2,5} base vertices
const SHAPES_3D = [
    ("ball",      Ball([0.0,0.0,0.0], 1.0),                                          :dodgerblue),
    ("cuboid",    Cuboid([0.0,0.0,0.0], [1.6,1.1,0.8], rotz(deg2rad(20))),           :tomato),
    ("ellipsoid", Ellipsoid([0.0,0.0,0.0], [1.3,0.9,0.6]),                           :mediumpurple),
    ("cylinder",  Cylinder([0.0,0.0,0.0], 0.8, 1.8, [0.25,0.1,1.0]),                 :seagreen),
    ("polyprism", PolygonalPrism(SVector(0.0,0.0,0.0), _pentagon, 1.6),              :darkorange),
    ("sectprism", SectoralPrism([0.0,0.0,0.0], 1.0, deg2rad(10), deg2rad(120), 1.6), :crimson),
]

@testset "visualization" begin
    @testset "2D shapes" begin
        @testset "$name" for (name, s, col) in SHAPES_2D
            fig = drawshape(s; color=(col, 0.55), axiskwargs=(; title="$name (2D)"))
            savefig("2d_$name.png", fig)
        end
    end

    @testset "3D shapes (perspective + slices)" begin
        @testset "$name" for (name, s, col) in SHAPES_3D
            lo, hi = bounds(s)
            ctr = (lo .+ hi) ./ 2

            # 3D perspective views from two angles.
            for (tag, az, el) in (("a", 1.275π, π/8), ("b", 0.30π, π/5))
                fig = drawshape(s; color=col, azimuth=az, elevation=el,
                                axiskwargs=(; title="$name (3D, view $tag)"))
                savefig("3d_$(name)_perspective_$tag.png", fig)
            end

            # 2D slices through the center along each Cartesian axis.
            for (ax, ci) in ((:x, ctr[1]), (:y, ctr[2]), (:z, ctr[3]))
                fig = drawshape(s(ax, ci); color=(col, 0.55),
                                axiskwargs=(; title="$name slice $ax=$(round(ci, digits=2))"))
                savefig("3d_$(name)_slice_$(ax).png", fig)
            end
        end
    end

    @testset "combined 2D scene" begin
        fig = Figure()
        ax = Axis(fig[1,1]; aspect=DataAspect(), title="combined 2D primitives")
        drawshape!(ax, Ball([-2.0,0.0], 1.0);                              color=(:dodgerblue,0.5))
        drawshape!(ax, Cuboid([1.0,0.5],[1.6,1.0], rot2(deg2rad(20)));     color=(:tomato,0.5))
        drawshape!(ax, Polygon{3}([0.0,-1.5], 1.0, deg2rad(90));           color=(:seagreen,0.5))
        drawshape!(ax, Sector([2.5,-1.5], 1.2, deg2rad(20), deg2rad(90));  color=(:darkorange,0.5))
        savefig("combined_2d.png", fig)
    end

    @testset "combined 3D scene" begin
        fig = Figure()
        ax = Axis3(fig[1,1]; aspect=:data, title="combined 3D primitives",
                   azimuth=1.275π, elevation=π/8)
        drawshape!(ax, Ball([-2.0,0.0,0.0], 0.9);                          color=:dodgerblue)
        drawshape!(ax, Cuboid([1.0,0.0,0.0],[1.2,1.0,0.8]);               color=:tomato)
        drawshape!(ax, Cylinder([0.0,2.0,0.0], 0.6, 1.6, [0.0,0.0,1.0]);  color=:seagreen)
        savefig("combined_3d.png", fig)
    end
end
