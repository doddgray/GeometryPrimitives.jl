# Makie visualization support for GeometryPrimitives, loaded automatically when Makie is
# available.  Provides `drawshape(s)` / `drawshape!(ax, s)` for every shape type:
#
#   * 2D shapes (`Shape2`, including `CrossSection`) are drawn as a filled boundary polygon
#     in an `Axis`.
#   * 3D shapes (`Shape3`) are drawn as a surface mesh in an `Axis3`.
#
# Both use the same idea: each primitive in this package is convex, hence star-shaped from
# any interior point, so its boundary in a given direction can be found by bisecting the
# level-set function `level` (≥ 0 inside, < 0 outside) along a ray from an interior point.
# This needs no external meshing dependency and renders with the pure-CPU CairoMakie
# backend.  For a 3D shape, a 2D slice is just `drawshape(s(:z, c))` (a `CrossSection`,
# which is itself a `Shape2`).
module GeometryPrimitivesMakieExt

import GeometryPrimitives as GP
using GeometryPrimitives: Shape2, Shape3, Cuboid, Prism, level, bounds
using Makie
using Makie: RGBAf, to_color
using Makie.GeometryBasics: GLTriangleFace
using LinearAlgebra: norm, dot, cross, normalize
using StaticArrays

# ---------------------------------------------------------------------------
# Boundary finding via bisection of the level-set along a ray from an interior point.
# ---------------------------------------------------------------------------

# Largest ρ ≥ 0 along p₀ + ρ*d̂ that is still inside s (level ≥ 0).  Assumes s is
# star-shaped from p₀ (true for the convex primitives here) so that the inside set along the
# ray is the interval [0, ρ].
function _boundary_radius(s, p₀::SVector{N}, d̂::SVector{N}, rmax::Real; iters::Int=60) where {N}
    # Make sure the far end is outside; expand if the initial guess was too small.
    hi = float(rmax)
    nexpand = 0
    while level(p₀ + hi*d̂, s) ≥ 0 && nexpand < 60
        hi *= 2; nexpand += 1
    end
    lo = zero(hi)  # level(p₀) ≥ 0 by assumption
    for _ in 1:iters
        mid = (lo + hi) / 2
        if level(p₀ + mid*d̂, s) ≥ 0
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

# Find an interior point (level > 0) of s within its bounding box.  Tries the box center
# first, then a coarse grid; returns `nothing` if none is found (e.g. an empty slice).
function _interior_point(s, lo::SVector{N}, hi::SVector{N}; n::Int=11) where {N}
    c = (lo + hi) / 2
    level(c, s) > 0 && return c
    best = c; bestlv = level(c, s)
    rng = ntuple(_ -> range(0, 1, length=n), Val(N))
    for idx in CartesianIndices(ntuple(_ -> n, Val(N)))
        t = SVector(ntuple(k -> rng[k][idx[k]], Val(N)))
        p = lo .+ t .* (hi - lo)
        lv = level(p, s)
        if lv > bestlv
            bestlv = lv; best = p
        end
    end
    return bestlv > 0 ? best : nothing
end

# ---------------------------------------------------------------------------
# 2D: filled boundary polygon.
# ---------------------------------------------------------------------------

# Ordered boundary points (counter-clockwise around an interior point) of a 2D shape.
function boundary_loop(s::Shape2; res::Integer=180)
    lo, hi = bounds(s)
    # Use the in-plane (first two) components: a genuine Shape2 already returns 2-vectors,
    # while a CrossSection reports its bounds in projected 3D coordinates (dims 1,2 are the
    # in-plane extent, dim 3 is along the slice normal).
    lo = SVector(lo[1], lo[2]); hi = SVector(hi[1], hi[2])
    p₀ = _interior_point(s, lo, hi)
    p₀ === nothing && return Point2f[]
    rmax = norm(hi - lo)
    pts = Vector{Point2f}(undef, res)
    for (i, θ) in enumerate(range(0, 2π, length=res+1)[1:res])
        d̂ = SVector(cos(θ), sin(θ))
        ρ = _boundary_radius(s, p₀, d̂, rmax)
        p = p₀ + ρ*d̂
        pts[i] = Point2f(p[1], p[2])
    end
    return pts
end

function GP.drawshape!(ax, s::Shape2; res::Integer=180, fill::Bool=true,
                       strokewidth::Real=1.5, kwargs...)
    pts = boundary_loop(s; res=res)
    isempty(pts) && return ax
    if fill
        poly!(ax, pts; strokewidth=strokewidth, kwargs...)
    else
        lines!(ax, push!(copy(pts), pts[1]); kwargs...)
    end
    return ax
end

function GP.drawshape(s::Shape2; res::Integer=180, axiskwargs=NamedTuple(), kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1]; aspect=DataAspect(), axiskwargs...)
    GP.drawshape!(ax, s; res=res, kwargs...)
    return fig
end

# ---------------------------------------------------------------------------
# 3D meshing and shading.
#
# CairoMakie applies no real lighting, so meshes are shaded manually: each surface point is
# tinted by a Lambert term n̂·L (ambient + diffuse) using a fixed light direction.  For the
# smooth shapes (Ball, Ellipsoid) the normal is approximated by the outward radial direction
# from the center; the polyhedral/extruded shapes (Cuboid, Prism and its Cylinder /
# PolygonalPrism / SectoralPrism aliases) are built as explicit flat-shaded faces with exact
# per-face normals.
# ---------------------------------------------------------------------------

const _LIGHT = normalize(SVector(0.4, 0.5, 0.85))
_lambert(n̂) = clamp(0.4 + 0.6*max(0.0, dot(n̂, _LIGHT)), 0.0, 1.0)

# Tint a base color toward black by the factor `t` ∈ [0,1], preserving its alpha.
function _shade(basecolor, t)
    c = to_color(basecolor)
    return RGBAf(t*c.r, t*c.g, t*c.b, c.alpha)
end

# Add a flat-shaded triangle soup `tris` (a list of (p₁,p₂,p₃)) to `ax`.  Each triangle is
# given a single color from the Lambert shading of its outward normal (oriented away from
# `center`), giving crisp faces.
function _flat_mesh!(ax, tris::Vector{NTuple{3,Point3f}}, center::Point3f, basecolor; kwargs...)
    nt = length(tris)
    verts = Vector{Point3f}(undef, 3nt)
    cols = Vector{RGBAf}(undef, 3nt)
    faces = Vector{GLTriangleFace}(undef, nt)
    for (t, (p1, p2, p3)) in enumerate(tris)
        n̂ = normalize(cross(p2 - p1, p3 - p1))
        dot(n̂, (p1 + p2 + p3)/3 - center) < 0 && (n̂ = -n̂)  # orient outward
        sh = _shade(basecolor, _lambert(n̂))
        i = 3t - 2
        verts[i], verts[i+1], verts[i+2] = p1, p2, p3
        cols[i] = cols[i+1] = cols[i+2] = sh
        faces[t] = GLTriangleFace(i, i+1, i+2)
    end
    mesh!(ax, verts, faces; color=cols, kwargs...)
    return ax
end

# Smooth ray-cast surface (UV-sphere topology) for star-shaped shapes; used for Ball and
# Ellipsoid (and as a generic fallback for any other Shape3).
function _uv_faces(nθ::Int, nφ::Int)
    faces = GLTriangleFace[]
    vidx(i, j) = (i-1)*nφ + ((j-1) % nφ) + 1
    for i in 1:nθ-1, j in 1:nφ
        a = vidx(i, j); b = vidx(i, j+1); c = vidx(i+1, j); d = vidx(i+1, j+1)
        push!(faces, GLTriangleFace(a, b, d))
        push!(faces, GLTriangleFace(a, d, c))
    end
    return faces
end

function surface_mesh(s::Shape3; nθ::Integer=64, nφ::Integer=96)
    lo, hi = bounds(s)
    lo = SVector{3}(lo); hi = SVector{3}(hi)
    all(isfinite, lo) && all(isfinite, hi) ||
        error("drawshape: shape has an infinite extent (e.g. a prism with infinite " *
              "height); construct it with a finite size before plotting.")
    p₀ = _interior_point(s, lo, hi)
    p₀ === nothing && error("drawshape: could not find an interior point of $s")
    rmax = norm(hi - lo)
    θs = range(0, π, length=nθ)
    φs = range(0, 2π, length=nφ+1)[1:nφ]
    verts = Vector{Point3f}(undef, nθ*nφ)
    k = 0
    for θ in θs, φ in φs
        d̂ = SVector(sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ))
        ρ = _boundary_radius(s, p₀, d̂, rmax)
        p = p₀ + ρ*d̂
        verts[k += 1] = Point3f(p[1], p[2], p[3])
    end
    return verts, _uv_faces(nθ, nφ)
end

function GP.drawshape!(ax, s::Shape3; color=:steelblue, nθ::Integer=64, nφ::Integer=96, kwargs...)
    verts, faces = surface_mesh(s; nθ=nθ, nφ=nφ)
    ctr = sum(verts) / length(verts)
    cols = [_shade(color, _lambert(normalize(v - ctr))) for v in verts]
    mesh!(ax, verts, faces; color=cols, kwargs...)
    return ax
end

# Exact flat-faced mesh for a 3D cuboid (8 corners → 6 quad faces → 12 triangles).
function _cuboid_tris(s::Cuboid{3})
    A = inv(s.p)  # columns are the (unit) cuboid axes
    corner(sx, sy, sz) = Point3f((s.c + A*SVector(sx*s.r[1], sy*s.r[2], sz*s.r[3]))...)
    p = [corner(sx, sy, sz) for sz in (-1,1) for sy in (-1,1) for sx in (-1,1)]
    idx(ix, iy, iz) = 1 + ix + 2iy + 4iz  # ix,iy,iz ∈ {0,1}
    quad(a, b, c, d) = ((p[a],p[b],p[c]), (p[a],p[c],p[d]))
    tris = NTuple{3,Point3f}[]
    append!(tris, quad(idx(0,0,0), idx(0,1,0), idx(0,1,1), idx(0,0,1)))  # x⁻
    append!(tris, quad(idx(1,0,0), idx(1,1,0), idx(1,1,1), idx(1,0,1)))  # x⁺
    append!(tris, quad(idx(0,0,0), idx(1,0,0), idx(1,0,1), idx(0,0,1)))  # y⁻
    append!(tris, quad(idx(0,1,0), idx(1,1,0), idx(1,1,1), idx(0,1,1)))  # y⁺
    append!(tris, quad(idx(0,0,0), idx(1,0,0), idx(1,1,0), idx(0,1,0)))  # z⁻
    append!(tris, quad(idx(0,0,1), idx(1,0,1), idx(1,1,1), idx(0,1,1)))  # z⁺
    return tris
end

function GP.drawshape!(ax, s::Cuboid{3}; color=:steelblue, kwargs...)
    return _flat_mesh!(ax, _cuboid_tris(s), Point3f(s.c...), color; kwargs...)
end

# Flat-faced mesh for a prism: extrude the 2D base boundary along the prism axis, with caps.
function _prism_tris(s::Prism; res::Integer=120)
    base = boundary_loop(s.b; res=res)  # base-plane boundary, ordered CCW
    length(base) < 3 && return NTuple{3,Point3f}[]
    Q = s.p'  # columns are the prism axes (u, v, and prism axis w)
    h = s.h2
    to3(uv, w) = Point3f((s.c + Q*SVector(uv[1], uv[2], w))...)
    ctr = sum(base) / length(base)  # 2D centroid (interior, since base is convex)
    tris = NTuple{3,Point3f}[]
    n = length(base)
    for i in 1:n
        j = i % n + 1
        bi, bj = base[i], base[j]
        a = to3(bi, -h); b = to3(bj, -h); c = to3(bj, h); d = to3(bi, h)
        push!(tris, (a, b, c)); push!(tris, (a, c, d))      # side wall
        push!(tris, (to3(ctr, -h), b, a))                   # bottom cap
        push!(tris, (to3(ctr,  h), d, c))                   # top cap
    end
    return tris
end

function GP.drawshape!(ax, s::Prism; color=:steelblue, res::Integer=120, kwargs...)
    return _flat_mesh!(ax, _prism_tris(s; res=res), Point3f(s.c...), color; kwargs...)
end

function GP.drawshape(s::Shape3; azimuth=1.275π, elevation=π/8,
                      axiskwargs=NamedTuple(), kwargs...)
    fig = Figure()
    ax = Axis3(fig[1,1]; aspect=:data, azimuth=azimuth, elevation=elevation, axiskwargs...)
    GP.drawshape!(ax, s; kwargs...)
    return fig
end

# ---------------------------------------------------------------------------
# `contour`/`heatmap` support for 2D shapes via the level-set grid (kept for
# compatibility): `contour(shp, (hres, vres))`.
# ---------------------------------------------------------------------------
function Makie.convert_arguments(P::GridBased, shp::Shape2, res::NTuple{2,Integer})
    lower, upper = bounds(shp)
    ∆ = upper - lower
    xs = range(lower[1] - GP.rtol(∆[1]), upper[1] + GP.rtol(∆[1]), length=res[1])
    ys = range(lower[2] - GP.rtol(∆[2]), upper[2] + GP.rtol(∆[2]), length=res[2])
    lvs = [level(SVector(x,y), shp) for x in xs, y in ys]
    return convert_arguments(P, xs, ys, lvs)
end

end # module
