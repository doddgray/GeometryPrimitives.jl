# Benchmarks for the primal evaluation and gradients of the differentiable API:
#   level(x, shape), surfpt_nearby(x, shape), volfrac(vxl, nout, r₀)
# with Enzyme.jl (forward and reverse mode) and Mooncake.jl (reverse mode), via
# DifferentiationInterface.jl, against central finite differences.
#
# Run with:  julia --project=benchmark benchmark/benchmarks.jl

using GeometryPrimitives
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using Random: MersenneTwister
using Printf

using DifferentiationInterface
using ADTypes: AutoEnzyme, AutoMooncake, AutoFiniteDifferences
import Enzyme, Mooncake, FiniteDifferences

const rng = MersenneTwister(7)

const BACKENDS = (
    "Enzyme reverse" => AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse),
                                     function_annotation=Enzyme.Const),
    "Enzyme forward" => AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward),
                                     function_annotation=Enzyme.Const),
    "Mooncake" => AutoMooncake(; config=nothing),
    "FiniteDifferences" => AutoFiniteDifferences(; fdm=FiniteDifferences.central_fdm(5, 1)),
)

rotmat2(θ) = SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ))
function circpts(K)
    ϕ = range(0.05, 1.95, length=K+1)
    return SMatrix{2,K}(ntuple(i -> isodd(i) ? cospi(ϕ[(i+1)÷2]) : sinpi(ϕ[i÷2]), Val(2K))...)
end

# (name, objective f(p::Vector)::Float64, p) triples covering representative shapes and the
# three differentiable functions, w.r.t. both query point x and shape parameters.
make_cases() = begin
    ball3 = Ball(SVector(0.1,0.2,0.3), 1.1)
    cub3 = Cuboid(SVector(0.1,0.2,0.3), SVector(1.0,1.2,1.4))
    pgon8 = Polygon(circpts(8))
    cyl = Cylinder(SVector(0.1,0.2,0.3), 0.8, 1.7, SVector(0.1,0.2,1.0))

    fvol(sh, N) = x -> begin
        xs = SVector{N}(x)
        δ = ones(SVector{N,eltype(x)})
        r₀, nout = surfpt_nearby(xs, sh)
        volfrac((xs - δ, xs + δ), nout, r₀)
    end

    return (
        ("level Ball3 / x",            x -> level(SVector{3}(x), ball3),                rand(rng, 3)),
        ("level Cuboid3 / x",          x -> level(SVector{3}(x), cub3),                 rand(rng, 3)),
        ("level Polygon8 / x",         x -> level(SVector{2}(x), pgon8),                rand(rng, 2)),
        ("level Cylinder / x",         x -> level(SVector{3}(x), cyl),                  rand(rng, 3)),
        ("surfpt Ball3 / x",           x -> sum(sum, surfpt_nearby(SVector{3}(x), ball3)), rand(rng, 3)),
        ("surfpt Cuboid3 / x",         x -> sum(sum, surfpt_nearby(SVector{3}(x), cub3)),  rand(rng, 3)),
        ("surfpt Polygon8 / x",        x -> sum(sum, surfpt_nearby(SVector{2}(x), pgon8)), rand(rng, 2)),
        ("surfpt Cylinder / x",        x -> sum(sum, surfpt_nearby(SVector{3}(x), cyl)),   rand(rng, 3)),
        ("volfrac Ball3 / x",          fvol(ball3, 3),                                  rand(rng, 3)),
        ("volfrac Cuboid3 / x",        fvol(cub3, 3),                                   rand(rng, 3)),
        ("level Ball3 / params",       p -> level(SVector(0.3,0.4,0.5), Ball(SVector(p[1],p[2],p[3]), p[4])),
                                       rand(rng, 4)),
        ("level Ellipsoid2 / params",  p -> level(SVector(0.3,0.4),
                                                  Ellipsoid(SVector(p[1],p[2]), SVector(p[3],p[4]), rotmat2(p[5]))),
                                       rand(rng, 5)),
        # Polygon5 (not a larger polygon): differentiating through the Polygon constructor
        # is the most expensive case to compile, especially with Enzyme reverse mode.
        ("surfpt Polygon5 / params",   p -> sum(sum, surfpt_nearby(SVector(0.3,0.4),
                                                  Polygon(0.1*SMatrix{2,5}(ntuple(k->p[k], Val(10))) + circpts(5)))),
                                       rand(rng, 10)),
        ("surfpt Cylinder / params",   p -> sum(sum, surfpt_nearby(SVector(0.3,0.4,0.5),
                                                  Cylinder(SVector(p[1],p[2],p[3]), p[4], p[5],
                                                           SVector(p[6],p[7],p[8])))),
                                       rand(rng, 8)),
    )
end

prettytime(t) = BenchmarkTools.prettytime(t)

function run_benchmarks(; seconds=1.0)
    cases = make_cases()
    results = []  # (case, primal_time, Dict(backend => (time, ratio)))
    for (name, f, p) in cases
        b_primal = @benchmark $f($p) seconds=seconds
        tp = minimum(b_primal).time
        row = Dict{String,Tuple{Float64,Float64}}()
        for (bname, backend) in BACKENDS
            prep = prepare_gradient(f, backend, zero(p))
            g = DifferentiationInterface.gradient(f, prep, backend, p)  # warm up + correctness sanity
            b = @benchmark DifferentiationInterface.gradient($f, $prep, $backend, $p) seconds=seconds
            tg = minimum(b).time
            row[bname] = (tg, tg/tp)
        end
        push!(results, (name, tp, row))
        @printf("%-28s primal %10s", name, prettytime(tp))
        for (bname, _) in BACKENDS
            tg, ratio = row[bname]
            @printf(" | %s %10s (%5.1fx)", bname, prettytime(tg), ratio)
        end
        println()
    end
    return results
end

function print_markdown(results)
    println("\n## Gradient benchmark results\n")
    header = "| case | primal |" * join((" $b |" for (b,_) in BACKENDS)) *
             "\n|---|---|" * repeat("---|", length(BACKENDS))
    println(header)
    for (name, tp, row) in results
        line = "| $name | $(prettytime(tp)) |"
        for (bname, _) in BACKENDS
            tg, ratio = row[bname]
            line *= @sprintf(" %s (%.1fx) |", prettytime(tg), ratio)
        end
        println(line)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run_benchmarks()
    print_markdown(results)
end
