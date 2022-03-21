using LinearAlgebra, StaticArrays, GeometryPrimitives, Test
circpts(n) = mapreduce(x->SVector{2}(reverse(sincospi(x)))',vcat,range(0.05,1.95,length=n+1)[1:n])
rand_verts(noise,n) = mapreduce(x->((1.0 + noise*rand())*SVector{2}(reverse(sincospi(x))))',vcat,range(0.05,1.95,length=n+1)[1:n])
using ChainRulesCore
ChainRulesCore.@non_differentiable circpts(::Any)
using ChainRules, Zygote, FiniteDifferences, ForwardDiff
# using CairoMakie
# xy1,v1 = SVector{2}(rand(2)),rand_verts(0.1,5)
# p1 = vcat(vec(xy1),vec(v1))
# lines(eachcol(circpts(20))...)
# lines(eachcol(SMatrix{10,2}(0.3*rand(20))+circpts(10))...)

function test_AD(f::Function,p;nFD=5)
    primal  =   f(p)
    gr_RM   =   first(Zygote.gradient(f,p))
    gr_FM   =   ForwardDiff.gradient(f,p)
    gr_FD   =   first(FiniteDifferences.grad(central_fdm(nFD,1),f,p))
    return isapprox(gr_RM,gr_FD,rtol=1e-4) && isapprox(gr_RM,gr_FD,rtol=1e-4)
end

function demo_shapes2D(p::Vector{T}=rand(17)) where T<:Real
    ε₁  =   diagm([p[1],p[1],p[1]])
    ε₂  =   diagm([p[2],p[3],p[3]])
    ε₃  =   diagm([1.0,1.0,1.0])
    b = Box(					# Instantiate N-D box, here N=2 (rectangle)
        p[4:5],					# c: center
        p[6:7],				# r: "radii" (half span of each axis)
        mapreduce(
            normalize,
            hcat,
            eachcol(reshape(p[8:11],(2,2))),
        ),			# axes: box axes
        ε₁,						# data: any type, data associated with box shape
        )

    s = Sphere(					# Instantiate N-D sphere, here N=2 (circle)
        p[12:13],					# c: center
        p[14],						# r: "radii" (half span of each axis)
        ε₂,						# data: any type, data associated with circle shape
        )

    t = regpoly(				# triangle::Polygon using regpoly factory method
        3,						# k: number of vertices
        p[15],					# r: distance from center to vertices
        π/2,					# θ: angle of first vertex
        p[16:17],					# c: center
        ε₃,						# data: any type, data associated with triangle
        )

	return ( t, s, b )
end

# 2D Sphere and Box functions
make_sphere2(p) = Sphere(SVector{2}(p[1:2]),p[3],p[4])
function make_box2(p)   # length(p) should be 9 
    return Box(					# Instantiate N-D box, here N=2 (rectangle)
        SVector{2}(p[1:2]),					# c: center
        SVector{2}(p[3:4]),				# r: "radii" (half span of each axis)
        SMatrix{2,2}(mapreduce(
            normalize,
            hcat,
            eachcol(reshape(p[5:8],(2,2)))),
        ),			# axes: box axes
        p[9],						# data: any type, data associated with box shape
    )
end
# 2D polygon functions
make_5gon(p) = Polygon(SMatrix{5,2}(reshape(p[1:10],(5,2)))*0.3+circpts(5),p[11])
make_10gon(p) = Polygon(SMatrix{10,2}(reshape(p[1:20],(10,2)))*0.1+circpts(10),p[21])
make_100gon(p) = Polygon(SMatrix{100,2}(reshape(p[1:200],(100,2)))*0.01+circpts(100),p[201])
make_8regpoly(p) = regpoly(Val(8),first(p),p[2],SVector{2}(p[3:4]),p[5])
make_regpoly(n,p) = regpoly(n,first(p),p[2],SVector{2}(p[3:4]),p[5])
fsh2D = ( make_sphere2, make_box2, make_5gon, make_10gon, make_8regpoly)
np_fsh2D = ( 4, 9, 11, 21, 5 )
# 3D Sphere and Box functions
# make_sphere3(p) = Sphere(SVector{3}(p[1:3]),p[4],p[5])
make_sphere3(p) = Sphere(SVector{3}(p[1],p[2],p[3]),p[4],p[5])
function make_box3(p)   # length(p) should be 9 
    return Box(					# Instantiate N-D box, here N=2 (rectangle)
        SVector{3}(p[1:3]),					# c: center
        SVector{3}(p[4:6]),				# r: "radii" (half span of each axis)
        SMatrix{3,3}(mapreduce(
            normalize,
            hcat,
            eachcol(reshape(p[7:15],(3,3)))),
        ),			# axes: box axes
        p[16],						# data: any type, data associated with box shape
    )
end
# 3D prism-related shape functions
make_6polyprism(p) = Prism(SVector{3}(p[14:16]),Polygon(SMatrix{6,2}(reshape(p[1:12],(6,2)))*0.3+circpts(6),p[13]),p[13])
make_cyl(p) = Cylinder(SVector{3}(p[1:3]),p[4],p[5],SVector{3}(p[6:8]),p[9])
fsh3D = ( make_sphere3, make_box3, make_6polyprism, make_cyl )
np_fsh3D = ( 5, 16, 16, 9 )
# make shapes with random parameters
# shapes2D = map((f,np)->f(rand(np)),zip(fsh2D,np_fsh2D))
# shapes3D = map((f,np)->f(rand(np)),zip(fsh3D,np_fsh3D))
# (sph2, bx2, pgn5, pgn10, pgn100, rp8) = shapes2D
# (sph3, bx3, ppr6, cylpr) = shapes3D

# # make 2D shapes
sph2 = make_sphere2(rand(4))
bx2 = make_box2(rand(9))
pgn5 = make_5gon(rand(11))
pgn10 = make_10gon(rand(21))
pgn100 = make_100gon(rand(201))
rp8 = make_8regpoly(rand(5))
rp60 = make_regpoly(60,rand(5))
shapes2D = ( sph2, bx2, pgn5, pgn10, rp8, )

# make 3D shapes
sph3 = make_sphere3(rand(5))
bx3 = make_box3(rand(16))
ppr6 = make_6polyprism(rand(16))
cylpr = make_cyl(rand(9))
shapes3D = ( sph3, bx3, ppr6, cylpr )
# shapes3D = map((f,np)->f(rand(np)),zip(fsh3D,np_fsh3D))
# (sph2, bx2, pgn5, pgn10, pgn100, rp8, rp60, sph3, bx3, ppr6, cylpr)map((f,np)->f(rand(np)),zip(vcat(fsh2D,fsh3D),vcat(np_fsh2D,np_fsh3D)))

@testset verbose = true "GeometryPrimitives AD Testing" begin
    @testset verbose = true "surfpt_nearby AD Gradients" begin
        @testset "2D x-vector gradients of surfpt_nearby(x,shape) for 2D shape type $(typeof(sh))" for sh in shapes2D
            p = rand(2)
            @test test_AD(x->sum(sum(surfpt_nearby(SVector{2}(x),sh))),p)
        end
        # @testset "3D x-vector gradients of surfpt_nearby(x,shape) for 2D shape type $(typeof(sh))" for sh in shapes2D
        #     p = rand(3)
        #     @test test_AD(x->sum(sum(surfpt_nearby(SVector{3}(x),sh))),p)
        # end
        @testset "3D x-vector gradients of surfpt_nearby(x,shape) for 3D shape type $(typeof(sh))" for sh in shapes3D
            p = rand(3)
            @test test_AD(x->sum(sum(surfpt_nearby(SVector{3}(x),sh))),p)
        end
        @testset "2D shape parameter gradients of surfpt_nearby(x::SVector{2},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
            xy = SVector{2}(rand(2))
            p = rand(np)
            @test test_AD(x->sum(sum(surfpt_nearby(xy,fsh(x)))),p)
        end
        # @testset "2D shape parameter gradients of surfpt_nearby(x::SVector{3},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
        #     xyz = SVector{3}(rand(3))
        #     p = rand(np)
        #     @test test_AD(x->sum(sum(surfpt_nearby(xyz,fsh(x)))),p)
        # end
        @testset "3D shape parameter gradients of surfpt_nearby(x::SVector{3},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh3D,np_fsh3D)
            xyz = SVector{3}(rand(3))
            p = rand(np)
            @test test_AD(x->sum(sum(surfpt_nearby(xyz,fsh(x)))),p)
        end
    end
    @testset verbose = true "volfrac AD Gradients" begin
        @testset "2D x-vector gradients of volfrac((x-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for 2D shape type $(typeof(sh))" for sh in shapes2D
            p = rand(2)
            @test test_AD(x->volfrac((SVector{2}(x)-SVector{2}(1.0,1.0),SVector{2}(x)+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(x),sh))...),p)
        end
        # @testset "3D x-vector gradients of volfrac((x-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for 2D shape type $(typeof(sh))" for sh in shapes2D
        #     p = rand(3)
        #     @test test_AD(x->volfrac((SVector{3}(x)-SVector{3}(1.0,1.0,1.0),SVector{3}(x)+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(x),sh))...),p)
        # end
        @testset "3D x-vector gradients of volfrac((x-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for 3D shape type $(typeof(sh))" for sh in shapes3D
            p = rand(3)
            @test test_AD(x->volfrac((SVector{3}(x)-SVector{3}(1.0,1.0,1.0),SVector{3}(x)+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(x),sh))...),p)
        end
        @testset "2D shape parameter gradients of volfrac((x::SVector{2}-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
            xy = SVector{2}(rand(2))
            p = rand(np)
            @test test_AD(x->volfrac((xy-SVector{2}(1.0,1.0),xy+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(xy,fsh(x)))...),p)
        end
        # @testset "2D shape parameter gradients of volfrac((x::SVector{3}-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
        #     xyz = SVector{3}(rand(3))
        #     p = rand(np)
        #     @test test_AD(x->volfrac((xyz-SVector{3}(1.0,1.0,1.0),xyz+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(xyz,fsh(x)))...),p)
        # end
        @testset "3D shape parameter gradients of volfrac((x::SVector{3}-δx,x+δx),reverse(surfpt_nearby(x,shape))...) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh3D,np_fsh3D)
            xyz = SVector{3}(rand(3))
            p = rand(np)
            @test test_AD(x->volfrac((xyz-SVector{3}(1.0,1.0,1.0),xyz+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(xyz,fsh(x)))...),p)
        end
    end
end

##
# Test Summary:                                                                                                                                                    | Pass  Total
# GeometryPrimitives AD Testing                                                                                                                                    |   36     36
#   surfpt_nearby AD Gradients                                                                                                                                     |   18     18
#   volfrac AD Gradients                                                                                                                                           |   18     18


