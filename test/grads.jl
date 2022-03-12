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
fsh2D = ( make_sphere2, make_box2, make_5gon, make_10gon, make_100gon, make_8regpoly)
np_fsh2D = ( 4, 9, 11, 21, 201, 5 )
# 3D Sphere and Box functions
make_sphere3(p) = Sphere(SVector{3}(p[1:3]),p[4],p[5])
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
make_6polyprism(p) = PolygonalPrism(SVector{3}(p[1:3]),SMatrix{6,2}(p[4:15]),p[16])
make_cyl(p) = Cylinder(SVector{3}(p[1:3]),p[4],p[5],SVector{3}(p[6:8]),p[9])
fsh3D = ( make_sphere3, make_box3, make_6polyprism, make_cyl )
np_fsh3D = ( 5, 16, 16, 9 )
# make shapes with random parameters
shapes2D = map((f,np)->f(rand(np)),zip(fsh2D,np_fsh2D))
shapes3D = map((f,np)->f(rand(np)),zip(fsh3D,np_fsh3D))
(sph2, bx2, pgn5, pgn10, pgn100, rp8) = shapes2D
(sph3, bx3, ppr6, cylpr) = shapes3D

# # make 2D shapes
# sph2 = make_sphere2(rand(4))
# bx2 = make_box2(rand(9))
# pgn5 = make_5gon(rand(11))
# pgn10 = make_10gon(rand(21))
# pgn100 = make_100gon(rand(201))
# rp8 = make_8regpoly(rand(5))
# rp60 = make_regpoly(60,rand(5))
# shapes2D = ( sph2, bx2, pgn5, pgn10, pgn100, rp8, rp60, )

# # make 3D shapes
# sph3 = make_sphere3(rand(5))
# bx3 = make_box3(rand(16))
# ppr6 = make_6polyprism(rand(16))
# cylpr = make_cyl(rand(9))
# shapes3D = ( sph3, bx3, ppr6, cylpr )
# shapes3D = map((f,np)->f(rand(np)),zip(fsh3D,np_fsh3D))
# (sph2, bx2, pgn5, pgn10, pgn100, rp8, rp60, sph3, bx3, ppr6, cylpr)map((f,np)->f(rand(np)),zip(vcat(fsh2D,fsh3D),vcat(np_fsh2D,np_fsh3D)))

Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{2}(p[12:13]),make_5gon(p[1:11])))),rand(13))
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[12:13]...,0.0)[1:2],make_5gon(p[1:11])))),rand(13))
test_AD(p->sum(sum(surfpt_nearby(SVector{2}(p[12:13]),make_5gon(p[1:11])))),rand(13))

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


##
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[6:8]),make_sphere3(p[1:5])))),rand(8))
test_AD(p->sum(sum(surfpt_nearby(SVector{3}(p[6:8]),make_sphere3(p[1:5])))),rand(8))
vf_sphere3(p) = volfrac((SVector{3}(p[6:8])-SVector{3}(0.1,0.1,0.1),SVector{3}(p[6:8])+SVector{3}(0.1,0.1,0.1)),reverse(surfpt_nearby(SVector{3}(p[6:8]),make_sphere3(p[1:5])))...)
Zygote.gradient(vf_sphere3,rand(8))
##
f = vf_sphere3
# p = [0., 0., 0., 1.0, 2.2,  0.99, 0.0, 0.0] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  1.01, 0.0, 0.0] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 0.99, 0.0] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 1.01, 0.0] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 0.0, 0.92] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 0.0, 1.01] + 0.01*rand(8)
p = [0., 0., 0., 1.0, 2.2,  0.0, 0.99/sqrt(2), 0.99/sqrt(2)] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 0.99/sqrt(2), 0.99/sqrt(2)] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.99/sqrt(2), 0.99/sqrt(2), 0.0 ] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.0, 0.99/sqrt(2), 0.99/sqrt(2)] + 0.01*rand(8)
# p = [0., 0., 0., 1.0, 2.2,  0.99/sqrt(3), 0.99/sqrt(3), 0.99/sqrt(3) ] + 0.01*rand(8)
nFD = 5
primal  =   f(p)
gr_RM   =   first(Zygote.gradient(f,p))
gr_FM   =   ForwardDiff.gradient(f,p)
gr_FD   =   first(FiniteDifferences.grad(central_fdm(nFD,1),f,p))
AD_test_ret = isapprox(gr_RM,gr_FD,rtol=1e-4) && isapprox(gr_RM,gr_FD,rtol=1e-4)
##
using GeometryPrimitives: corner_bits, rvol_quadsect, rvol_gensect, isquadsect, edgedir_quadsect
s = make_sphere3(p[1:5])
x = SVector{3}(p[6:8])
vxl = (x-SVector{3}(0.1,0.1,0.1),x+SVector{3}(0.1,0.1,0.1))
r₀, nout = surfpt_nearby(x,s)
nr₀ = nout⋅r₀
cbits, n_on = corner_bits(vxl, nout, nr₀)
n_in = count_ones(cbits)
isquadsect(cbits)
w = edgedir_quadsect(cbits)
if isequal(w,1)
    u, v, _w = (2, 3, 1)
elseif isequal(w,2)
    u, v, _w = (3, 1, 2)
elseif isequal(w,3)
    u, v, _w = (1,2,3)
else
    println("bad return value from edgedir_quadsect in rvol_quadsect: w = $w")
end
∆w = vxl[2][w] - vxl[1][w]
nu, nv, nw = nout[u], nout[v], nout[w]
mean_cepts1 = 4*nr₀
for sv in (1,2), su in (1,2)
    mean_cepts1 -= nu*vxl[su][u] + nv*vxl[sv][v]
end
mean_cepts1 /=  nw * 4*∆w
mean_cepts = ( 4*nr₀ - ( nu*vxl[1][u] + nv*vxl[1][v] ) - ( nu*vxl[2][u] + nv*vxl[1][v] ) - ( nu*vxl[1][u] + nv*vxl[2][v] ) - ( nu*vxl[2][u] + nv*vxl[2][v] ) ) / ( nw * 4*∆w ) 
sw = nw>0 ? 1 : 2
rvqs1 = abs(mean_cepts - vxl[sw][w]/∆w)
rvolq = rvol_quadsect(vxl, nout, nr₀, cbits)
rvolg = rvol_gensect(vxl, nout, nr₀, cbits)

##
X, Y, Z = 1, 2, 3
XYZ, YZX, ZXY = (X,Y,Z), (Y,Z,X), (Z,X,Y)
UVW = YZX, ZXY, XYZ
N, P = 1, 2  # negative, positive
NP = (N, P)


##

Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[17:19]),make_box3(p[1:16])))),rand(19))
test_AD(p->sum(sum(surfpt_nearby(SVector{3}(p[17:19]),make_box3(p[1:16])))),rand(19))
Zygote.gradient(p->volfrac((SVector{3}(p[17:19])-SVector{3}(1.0,1.0,1.0),SVector{3}(p[17:19])+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(p[17:19]),make_box3(p[1:16])))...),rand(19))
test_AD(p->volfrac((SVector{3}(p[17:19])-SVector{3}(1.0,1.0,1.0),SVector{3}(p[17:19])+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(p[17:19]),make_box3(p[1:16])))...),rand(19))



##
p0 = rand(20)

sh1 = demo_shapes2D(p0);
pgn1,sph1,bx1 = sh1;
xy1 = SVector{2}(p0[18:19])
xyz1 = SVector{3}(p0[18:20])

@test isnothing(Zygote.gradient(in,xy1,pgn1))
@test isnothing(Zygote.gradient(in,xyz1,pgn1))
@test isa(surfpt_nearby(xy1,pgn1),Tuple{SVector{2, Float64}, SVector{2, Float64}})
@test isa(surfpt_nearby(xyz1,pgn1),Tuple{SVector{2, Float64}, SVector{2, Float64}})


# Test functions for differentiability of GeometryPrimitives methods for various Shape types and 2D and 3D vector inputs
f_spn_21(p) = sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[1])))
f_spn_22(p) = sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[2])))
f_spn_23(p) = sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3])))

f_spn_31(p) = sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[1])))
f_spn_32(p) = sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[2])))
f_spn_33(p) = sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[3])))

f_volfrac21(p) = volfrac((SVector{2}(p[18:19])-SVector{2}(1.0,1.0),SVector{2}(p[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[1]))...)
f_volfrac22(p) = volfrac((SVector{2}(p[18:19])-SVector{2}(1.0,1.0),SVector{2}(p[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[2]))...)
f_volfrac23(p) = volfrac((SVector{2}(p[18:19])-SVector{2}(1.0,1.0),SVector{2}(p[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3]))...)

f_volfrac31(p) = volfrac((SVector{3}(p[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[1]))...)
f_volfrac32(p) = volfrac((SVector{3}(p[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[2]))...)
f_volfrac33(p) = volfrac((SVector{3}(p[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[3]))...)

test_fns = [ 
    f_spn_21, f_spn_22, f_spn_23, 
    f_spn_31, f_spn_32, f_spn_33, 
    f_volfrac21, f_volfrac22, f_volfrac23, 
    f_volfrac31, f_volfrac32, f_volfrac33, 
]



test_AD(p0) do p
    sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3])))
end

sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3])))
p0 = rand(20)
p = rand(20)
demo_shapes2D(p)[1] |> typeof
demo_shapes2D(p)[2] |> typeof
demo_shapes2D(p)[3] |> typeof
sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[1])))
sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[2])))
sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3])))


sort_v_if_needed(SMatrix{10,2}(p10gon[1:20]))
sort_v_if_needed(SMatrix{10,2}(rand(20)))

p10gon = [sort_v_if_needed(SMatrix{10,2}(rand(20)))...,rand(4)...]
# p10g1 = Polygon(sort_v_if_needed(SMatrix{10,2}(p10gon[1:20])),p10gon[21])
circpts(n) = mapreduce(x->SVector{2}(reverse(sincospi(x)))',vcat,range(0,-2,length=n+1)[1:n])
# lines(eachcol(circpts(20))...)
lines(eachcol(SMatrix{10,2}(0.3*rand(20))+circpts(10))...)
p10g1 = Polygon(SMatrix{10,2}(0.3*rand(20))+circpts(10),rand())
cp10 = circpts(10)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{2}(p[22:23]),Polygon(SMatrix{10,2}(0.3*p[1:20])+cp10,p[21])))),rand(23))
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[18:20]),make_5gon(p)))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[1]))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[1]))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[2]))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[2]))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{2}(p[18:19]),demo_shapes2D(p)[3]))),p0)
Zygote.gradient(p->sum(sum(surfpt_nearby(SVector{3}(p[18:20]),demo_shapes2D(p)[3]))),p0)

volfrac((SVector{2}(p0[18:19])-SVector{2}(1.0,1.0),SVector{2}(p0[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p0[18:19]),demo_shapes2D(p0)[1]))...)
volfrac((SVector{2}(p0[18:19])-SVector{2}(1.0,1.0),SVector{2}(p0[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p0[18:19]),demo_shapes2D(p0)[2]))...)
volfrac((SVector{2}(p0[18:19])-SVector{2}(1.0,1.0),SVector{2}(p0[18:19])+SVector{2}(1.0,1.0)),reverse(surfpt_nearby(SVector{2}(p0[18:19]),demo_shapes2D(p0)[3]))...)

volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[1]))...)
volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[2]))...)
volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,0.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,0.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[3]))...)

volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,1.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[1]))...)
volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,1.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[2]))...)
volfrac((SVector{3}(p0[18:20])-SVector{3}(1.0,1.0,1.0),SVector{3}(p0[18:20])+SVector{3}(1.0,1.0,1.0)),reverse(surfpt_nearby(SVector{3}(p0[18:20]),demo_shapes2D(p0)[3]))...)



gr_volfrac1_RM = Zygote.gradient(f_volfrac1,p0)[1]
gr_volfrac1_FD = FiniteDifferences.grad(central_fdm(5,1),f_volfrac1,p0)[1]
@show isapprox(gr_volfrac1_RM,gr_volfrac1_FD,rtol=1e-4)


gr_volfrac2_RM = Zygote.gradient(f_volfrac2,p0)[1]
gr_volfrac2_FD = FiniteDifferences.grad(central_fdm(5,1),f_volfrac2,p0)[1]
@show isapprox(gr_volfrac2_RM,gr_volfrac2_FD,rtol=1e-4)


gr_volfrac3_RM = Zygote.gradient(f_volfrac3,p0)[1]
gr_volfrac3_FD = FiniteDifferences.grad(central_fdm(5,1),f_volfrac3,p0)[1]
@show isapprox(gr_volfrac3_RM,gr_volfrac3_FD,rtol=1e-4)

## prototype fixes to Box "max deviation" `m` vector calculation that breaks Zygote AD
p_inv2=SMatrix{2,2}(rand(2,2))
r2=SVector{2}(rand(2))
p_inv3=SMatrix{3,3}(rand(3,3))
r3=SVector{3}(rand(3))
using GeometryPrimitives: signmatrix
@show m2_old = let p_inv=p_inv2, r=r2
    smr = signmatrix(r)
    A = p_inv .* r'
    m = maximum(abs.(A * smr), dims=2)[:,1]
end
@show m2_new1 = let p_inv=p_inv2, r=r2
    smr = signmatrix(r)
    @tullio (max) m[i] := p_inv[i,k] * r[i] * smr[i,j]
end
@show m2_new2 = let p_inv=p_inv2, r=r2
    smr = signmatrix(r)
    A = mapreduce(*,hcat,eachcol(p_inv),r)
    m_in = map(abs,A*smr)
    @tullio (max) m_out[i] := m_in[i,j]
end
@show m2_new3 = let p_inv=p_inv2, r=r2
    smr = signmatrix(r)
    A = p_inv .* r'
    m_in = abs.(A * smr)
    @tullio (max) m_out[i] := m_in[i,j]
end

@show m3_old = let p_inv=p_inv3, r=r3
    smr = signmatrix(r)
    A = p_inv .* r'
    m = maximum(abs.(A * smr), dims=2)[:,1]
end
@show m3_new1 = let p_inv=p_inv3, r=r3
    smr = signmatrix(r)
    @tullio (max) m[i] := p_inv[i,k] * r[i] * smr[i,j]
end
@show m3_new2 = let p_inv=p_inv3, r=r3
    smr = signmatrix(r)
    A = mapreduce(*,hcat,eachcol(p_inv),r)
    m_in = map(abs,A*smr)
    @tullio (max) m_out[i] := m_in[i,j]
end
@show m3_new3 = let p_inv=p_inv3, r=r3
    smr = signmatrix(r)
    A = p_inv .* r'
    m_in = abs.(A * smr)
    @tullio (max) m_out[i] := m_in[i,j]
end
## prototype fixes to regpoly constructor vertex matrix `v` calculation that breaks Zygote AD
K = rand(5:100)
r = rand()  # distance between center and each vertex
θ = rand()  # angle from +x-direction towards first vertex; π/2 corresponds to +y-direction
c = SVector{2}(rand(2)) # center location
SVector(ntuple(k->k-1, Val(K)))
@show v_old = let K=K,r=r,θ=θ,c=c
    ∆θ = 2π / K
    θs = θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    v = c' .+ r .* [cos.(θs) sin.(θs)]  # SMatrix{K,2}: locations of vertices
end

@show v_new = let K=K,r=r,θ=θ,c=c
    ∆θ = 2π / K
    # θs = θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    θs = (θ * ones(K)) + collect(0:(K-1))*2π/K  
    # v = c' + (r * [cos.(θs) sin.(θs)])  # SMatrix{K,2}: locations of vertices
    v = mapreduce(vcat,θs) do tt
        SMatrix{1,2}(reverse(sincos(tt)))*r + c'
    end
end
@show v_new2 = let K=K,r=r,θ=θ,c=c
    ∆θ = 2π / K
    # θs = θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    θs = (θ * ones(K)) + collect(0:(K-1))*2π/K  
    # v = c' + (r * [cos.(θs) sin.(θs)])  # SMatrix{K,2}: locations of vertices
    v = reduce(vcat,map(tt->SMatrix{1,2}(reverse(sincos(tt)))*r + c', θs)) 
end
@assert v_old ≈ v_new
@assert v_old ≈ v_new2
## prototype fixes to mutating `sort_v_if_needed` function in Polygon constructor that breaks Zygote AD
K = rand(5:20)
v_in = SMatrix{K,2}(rand(2K))
using Statistics: mean
using GeometryPrimitives: sort_v_if_needed, ∆ϕ, n_norm
v = sort_v_if_needed(v_in)
∆v = v - circshift(v,1)

lines(vcat(v[:,1],v[1,1]),vcat(v[:,2],v[1,2]))
convex_inds1 = findall(x->x>0,∆ϕ(v-circshift(SMatrix{K,2}(v),1)))
v1 = v[convex_inds1,:]
convex_inds2 = findall(x->x>0,∆ϕ(v1-circshift(SMatrix{49,2}(v1),1)))
v2 = v1[convex_inds2,:]
convex_inds3 = findall(x->x>0,∆ϕ(v2-circshift(SMatrix{K,2}(v2),1)))
v3 = v1[convex_inds3,:]
lines(vcat(v[convex_inds,1],v[first(convex_inds),1]),vcat(v[convex_inds,2],v[first(convex_inds),2]))
all(∆ϕ(∆v) .> 0)
v = copy(v_in)

w1 = v .- mean(v,dims=1)
vm = sum(v,dims=1)/K
w2 = mapreduce(x->x'-vm,vcat,eachrow(v))
ϕ1 = mod.(atan.(w1[:,2], w1[:,1]), 2π)
ϕ2 = map(x->mod(atan(last(x),first(x)),2π),eachrow(w2))
@assert w1 ≈ w2
@assert ϕ2 ≈ ϕ1

Δv = v - circshift(v,1)
n_older = [Δv[:,2] -Δv[:,1]] ./ norm.(eachrow(Δv))
n_old = n_norm(Δv)
M_signs = SMatrix{2,2}(0,-1, 1,0)
n_new = mapreduce(x->normalize(M_signs*x)',vcat,eachrow(Δv))
@assert n_older ≈ n_old
@assert n_older ≈ n_new
@assert n_old ≈ n_new
@assert v_old ≈ v_new2
## Polygon constructor still breaking Zygote, must debug internals
using LinearAlgebra, StaticArrays, GeometryPrimitives, Test
using ChainRulesCore
using ChainRules, Zygote, FiniteDifferences, ForwardDiff
using GeometryPrimitives: polygon_vert_sortperms, l_bnds_poly, u_bnds_poly, ∆ϕ, n_norm
using Tullio
using CairoMakie
function test_AD(f::Function,p;nFD=5)
    primal  =   f(p)
    gr_RM   =   first(Zygote.gradient(f,p))
    gr_FM   =   ForwardDiff.gradient(f,p)
    gr_FD   =   first(FiniteDifferences.grad(central_fdm(nFD,1),f,p))
    return isapprox(gr_RM,gr_FD,rtol=1e-4) && isapprox(gr_RM,gr_FD,rtol=1e-4)
end
circpts(n) = mapreduce(x->SVector{2}(reverse(sincospi(x)))',vcat,range(0.05,1.95,length=n+1)[1:n])
rand_verts(noise,n) = mapreduce(x->((1.0 + noise*rand())*SVector{2}(reverse(sincospi(x))))',vcat,range(0.05,1.95,length=n+1)[1:n])
xy1,v1 = SVector{2}(rand(2)),rand_verts(0.1,5)
p1 = vcat(vec(xy1),vec(v1))

function fpgn1(xy::SVector{2,T1},v_in::SMatrix{K,2,T2}) where {K,T1<:Real,T2<:Real}
    # v = sort_v_if_needed(v_in)
    clockwise_sortperms = polygon_vert_sortperms(v_in)
    if isequal(clockwise_sortperms,collect(1:K))
        v = v_in
    else
        v = v_in[clockwise_sortperms,:]
	end
    ∆v = v - circshift(v,1)
    @assert all(∆ϕ(∆v) .> 0) # "v = $v must represent vertices of convex polygon, but v seems non-convex"
    n0 = ∆v * [	 0.     -1.     ;   1.     0.  ] # = [∆v[:,2] -∆v[:,1]]  # outward normal directions to edges
	@tullio nnorm[k] := ∆v[k,j]^2 |> sqrt
	@tullio n[k,j] := n0[k,j] / nnorm[k] # normalize
	l = l_bnds_poly(v)
	u = u_bnds_poly(v)
	sz = u-l
    rbnd = max(sz.data[1],sz.data[2])*Base.rtoldefault(T2)
    pgn = Polygon{K,2K,Float64,T2}(v,SMatrix{K,2,T2}(n),l,u,sz,rbnd,1.1) # Polygon{K,2K,D}(v,n,data)
    return sum(sum(surfpt_nearby(xy,pgn)))
end

xy1,v1 = SVector{2}(rand(2)),rand_verts(0.1,5)
p1 = vcat(vec(xy1),vec(v1))
res1 = fpgn1(xy1,v1)
gr1 = Zygote.gradient(fpgn1,xy1,v1)
test_AD(p->fpgn1(SVector{2}(p[1:2]),SMatrix{5,2}(reshape(p[3:12],(5,2)))), p1)

function fpgn2(xy::SVector{2,T1},v::SMatrix{K,2,T2}) where {K,T1<:Real,T2<:Real}
    ∆v = v - circshift(v,1)
    @assert all(∆ϕ(∆v) .> 0) # "v = $v must represent vertices of convex polygon, but v seems non-convex"
    n0 = ∆v * [	 0.     -1.     ;   1.     0.  ] # = [∆v[:,2] -∆v[:,1]]  # outward normal directions to edges
	@tullio nnorm[k] := ∆v[k,j]^2 |> sqrt
	@tullio n[k,j] := n0[k,j] / nnorm[k] # normalize
	l = l_bnds_poly(v)
	u = u_bnds_poly(v)
	sz = u-l
    rbnd = max(sz.data[1],sz.data[2])*Base.rtoldefault(T2)
    pgn = Polygon{K,2K,Float64,T2}(v,SMatrix{K,2,T2}(n),l,u,sz,rbnd,1.1) # Polygon{K,2K,D}(v,n,data)
    return volfrac((xy-[1.0,1.0],xy+[1.0,1.0]),reverse(surfpt_nearby(xy,pgn))...)
end
res2 = fpgn2(xy1,v1)
gr2 = Zygote.gradient(fpgn2,xy1,v1)
test_AD(p->fpgn2(SVector{2}(p[1:2]),SMatrix{5,2}(reshape(p[3:12],(5,2)))), p1)

## Polygon AD debug scratch
v1 = rand_verts(0.3,20)
clockwise_perm = sortperm(collect(eachrow(v1)),by=x->mod(atan(last(x),first(x)),2π))
isequal(clockwise_perm,collect(1:size(v1,1)))
v1sort = sort_v_if_needed(copy(v1))
fig,ax,lv1 = lines(eachcol(v1)...)
ax2 = fig[1,2] = Axis(fig)
lv2 = lines!(ax2,eachcol(v1sort)...); fig
# maybe better ways to sort vertices for clockwise order
# sortslices(v1,dims=1,by=x->mod(atan(last(x),first(x)),2π))
dphi1 = ∆ϕ((v1-circshift(v1,1)))
∆ϕ((v1-circshift(v1,1))) ≈  Δϕ1

w = v1 .- mean(v1, dims=1)  # v in center-of-mass coordinates
ϕ = mod.(atan.(w[:,2], w[:,1]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
issorted(ϕ)
v0 = circpts(10)
w0 = v0 .- mean(v0, dims=1)  # v in center-of-mass coordinates
ϕ0 = mod.(atan.(w0[:,2], w0[:,1]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
issorted(ϕ0)

K=size(v1,1)
Δv1 = (v1-circshift(v1,1))
Δz = Δv1[:,1] + im * Δv1[:,2]  # SVector{K}: edge directions as complex numbers
icurr = ntuple(identity, Val(K-1))
inext = ntuple(x->x+1, Val(K-1))
Δϕ1 =  angle.(Δz[SVector(inext)] ./ Δz[SVector(icurr)]) 
Δϕ11 =  map(x->atan(last(x),first(x)), eachrow(v1) )
Δϕ12 = circshift(Δϕ11,1) - Δϕ11
pgn1 = Polygon(v1,rand())
Δϕ2 = angle.(Δz)
zip(Δz,circshift(Δz,1)) |> first |> zpair->angle(last(zpair)/first(zpair))
Δϕ3 = map(zpair->angle(first(zpair)/last(zpair)),zip(Δz,circshift(Δz,1)))[2:end]
Δϕ4 = map(zpair->angle(first(zpair)/last(zpair)),zip(Δz[2:20],Δz[1:19]))
vz1 = SVector{2,ComplexF64}(1.0,1.0im)
Δϕ5 = map(xypair->angle((first(xypair)*vz1 )/(last(xypair)*vz1)),zip(eachrow(Δv1)[2:20],eachrow(Δv1)[1:19]))
@assert Δϕ1 ≈ Δϕ3
@assert Δϕ1 ≈ Δϕ4
@assert sort_v_if_needed(copy(v1)) ≈ v1

##

function ff1(p)
    shapes = demo_shapes2D(p)
    data_sum = sum( map(ss->getproperty(ss,:data),shapes) ) 
    return sum(data_sum)
end

@show valff1 = ff1(p0)
@show grff1 = Zygote.gradient(ff1,p0)[1]
@show typeof(grff1)

# ff2 series, variables shapes and constant points

function ff212(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(0.1,0.2)
    return first(first(surfpt_nearby(xyz,shapes[1])))
end
@show valff212 = ff212(p0)
@show grff212 = Zygote.gradient(ff212,p0)[1]
@show typeof(grff212)

function ff213(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(0.1,0.2,0.3)
    return first(first(surfpt_nearby(xyz,shapes[1])))
end
@show valff213 = ff213(p0)
@show grff213 = Zygote.gradient(ff213,p0)[1]
@show typeof(grff213)

function ff222(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(0.1,0.2)
    return first(first(surfpt_nearby(xyz,shapes[2])))
end
@show valff222 = ff222(p0)
@show grff222 = Zygote.gradient(ff222,p0)[1]
@show typeof(grff222)

function ff223(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(0.1,0.2,0.3)
    return first(first(surfpt_nearby(xyz,shapes[2])))
end
@show valff223 = ff223(p0)
@show grff223 = Zygote.gradient(ff223,p0)[1]
@show typeof(grff223)

function ff232(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(0.1,0.2)
    return first(first(surfpt_nearby(xyz,shapes[3])))
end
@show valff232 = ff232(p0)
@show grff232 = Zygote.gradient(ff232,p0)[1]
@show typeof(grff232)

function ff233(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(0.1,0.2,0.3)
    return first(first(surfpt_nearby(xyz,shapes[3])))
end
@show valff233 = ff233(p0)
@show grff233 = Zygote.gradient(ff233,p0)[1]
@show typeof(grff233)

# ff3 series, now with variable shapes and points

p0 = rand(20)

function ff312(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(p[18],p[19])
    return first(first(surfpt_nearby(xyz,shapes[1])))
end
@show valff312 = ff312(p0)
@show grff312 = Zygote.gradient(ff312,p0)[1]
@show typeof(grff312)

function ff313(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(p[18],p[19],p[20])
    return first(first(surfpt_nearby(xyz,shapes[1])))
end
@show valff313 = ff313(p0)
@show grff313 = Zygote.gradient(ff313,p0)[1]
@show typeof(grff313)

function ff322(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(p[18],p[19])
    return first(first(surfpt_nearby(xyz,shapes[2])))
end
@show valff322 = ff322(p0)
@show grff322 = Zygote.gradient(ff322,p0)[1]
@show typeof(grff322)

function ff323(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(p[18],p[19],p[20])
    return first(first(surfpt_nearby(xyz,shapes[2])))
end
@show valff323 = ff323(p0)
@show grff323 = Zygote.gradient(ff323,p0)[1]
@show typeof(grff323)

function ff332(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{2,Float64}(p[18],p[19])
    return first(first(surfpt_nearby(xyz,shapes[3])))
end
@show valff332 = ff332(p0)
@show grff332 = Zygote.gradient(ff332,p0)[1]
@show typeof(grff332)

function ff333(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(p[18],p[19],p[20])
    return first(first(surfpt_nearby(xyz,shapes[3])))
end
@show valff333 = ff333(p0)
@show grff333 = Zygote.gradient(ff333,p0)[1]
@show typeof(grff333)

##


function ff2(p)
    shapes = demo_shapes2D(p)
    xyz = SVector{3,Float64}(0.0,0.1,0.2)
    data_sum = sum( map(ss->first(first(surfpt_nearby(xyz,ss))),shapes) ) 
    return data_sum
end
ff2(pff1)
grff2 = Zygote.gradient(ff2,pff1)[1]
typeof(grff2)


shapes = demo_shapes2D(rand(17))
ff31(x,y,z) = first(first(surfpt_nearby(SVector{3}(x,y,z),shapes[1])))
ff31(0.3,0.2,1.1)
Zygote.gradient(ff31,0.3,0.2,1.1)

ff312(x,y,z) = first(first(surfpt_nearby(SVector{2}(x,y),shapes[1])))
ff312(0.3,0.2,1.1)
Zygote.gradient(ff312,0.3,0.2,1.1)

ff322(x,y,z) = first(first(surfpt_nearby(SVector{2}(x,y),3])))
ff322(0.3,0.2,1.1)
Zygote.gradient(ff322,0.3,0.2,1.1)

ff32(x,y,z) = first(first(surfpt_nearby(SVector{3}(x,y,z),shapes[2])))
ff32(0.3,0.2,1.1)
Zygote.gradient(ff32,0.3,0.2,1.1)















##

function foo2(p)
    shapes = demo_shapes2D(p)
    matvals = (getproperty.(shapes, (:data,))..., diagm([0.9,0.9,0.9]),)
    crnrs = (SVector{2,Float64}(0.1,0.1), SVector{2,Float64}(0.2,0.2), SVector{2,Float64}(0.1,0.25), SVector{2,Float64}(-0.1,0.14) )
    return smoov1_single(shapes,matvals,minds,crnrs)
end

foo2s(p) = sum(foo2(p))

function foo3(p)
    shapes = demo_shapes2D(p)
    matvals = (getproperty.(shapes, (:data,))..., diagm([0.9,0.9,0.9]),)
    crnrs = (SVector{3,Float64}(0.1,0.1,0.1), SVector{3,Float64}(0.2,0.2,0.2), SVector{3,Float64}(0.1,0.3,0.25), SVector{3,Float64}(-0.1,0.13,0.14) )
    return smoov1_single(shapes,matvals,minds,crnrs)
end

foo3s(p) = sum(foo3(p))

p0 = rand(17)
foo2(p0)
foo2s(p0)
Zygote.gradient(foo2s,p0)

grsm1s1          =   Zygote.gradient(p2) do p
    shapes2 = demo_shapes2D(p)
    matvals2 = (getproperty.(shapes2, (:data,))..., diagm([0.9,0.9,0.9]),)
    crnrs = (SVector{3,Float64}(0.1,0.1,0.1), SVector{3,Float64}(0.2,0.2,0.2), SVector{3,Float64}(0.1,0.3,0.25), SVector{3,Float64}(-0.1,0.13,0.14) )
    smval = smoov1_single(shapes2,matvals2,minds2,crnrs)
    return sum(smval)
end

grsm1s2          =   Zygote.gradient(p2) do p
    grid  =   Grid(3.0,4.0,64,128)
    shapes2 = demo_shapes2D(p)
    matvals2 = (getproperty.(shapes2, (:data,))..., diagm([0.9,0.9,0.9]),)
    smval = smoov1_single(shapes2,matvals2,minds2,first(corners(grid)))
    return sum(smval)
end




##
using Tullio
using GeometryPrimitives: signmatrix
p = rand(20)
c = SVector{3}(p[1:3])					# c: center
d = SVector{3}(p[4:6])				# r: "radii" (half span of each axis)
axes = SMatrix{3,3}(mapreduce(
    normalize,
    hcat,
    eachcol(reshape(p[7:15],(3,3)))),
)			# axes: box axes
data = p[16]		
r = 0.5d
@tullio axnorm[i] := axes[i,j]^2 |> sqrt
@tullio p_inv[i,j] := axes[i,j] / axnorm[j]
# p_inv = axes ./ sqrt.(sum(abs2,axes,dims=1))
smr = signmatrix(r)

A = p_inv .* r'
@show m1 = maximum(abs.(A * smr), dims=2)[:,1]
@show @tullio (max) m[j] := p_inv[i,k] * r[i] * smr[i,j] nograd=smr
A2 = [-1.1  2.2 ; -3.3  4.4]
abs(A2)
@tullio (max) m[j] := p_inv[i,k] * r[i] * smr[i,j] nograd=smr
bx_out = Box{3,9,Float64,Float64}(c, 0.5d, SMatrix(inv(p_inv)),c-SVector(m),c+SVector(m),data)
