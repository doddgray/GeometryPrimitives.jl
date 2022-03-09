using LinearAlgebra, StaticArrays, GeometryPrimitives, Test
using ChainRules, Zygote, FiniteDifferences, ForwardDiff
using GeometryPrimitives: sort_v_if_needed

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
make_5gon(p) = Polygon(sort_v_if_needed(SMatrix{5,2}(p[1:10])),p[11])
make_10gon(p) = Polygon(sort_v_if_needed(SMatrix{10,2}(p[1:20])),p[21])
make_100gon(p) = Polygon(sort_v_if_needed(SMatrix{100,2}(p[1:200])),p[201])
make_8regpoly(p) = regpoly(Val(8),first(p),p[2],SVector{2}(p[3:4]),p[5])
make_regpoly(n,p) = regpoly(n,first(p),p[2],SVector{2}(p[3:4]),p[5])
fsh2D = ( make_sphere2, make_box2, make_5gon, make_10gon, make_100gon, make_8regpoly)
np_fsh2D = ( 11, 21, 201, 5, 5 )
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
make_6polyprism(p) = PolygonalPrism(SVector{3}(p[1:3]),sort_v_if_needed(SMatrix{6,2}(p[4:15])),p[16])
make_cyl(p) = Cylinder(SVector{3}(p[1:3]),p[4],p[5],SVector{3}(p[6:8]),p[9])
fsh3D = ( make_sphere3, make_box3, make_6polyprism, make_cyl )
np_fsh3D = ( 5, 16, 16, 9 )
# make 2D shapes
sph2 = make_sphere2(rand(4))
bx2 = make_box2(rand(9))
pgn5 = make_5gon(rand(11))
pgn10 = make_10gon(rand(21))
pgn100 = make_100gon(rand(201))
rp8 = make_8regpoly(rand(5))
rp60 = make_regpoly(60,rand(5))
shapes2D = ( sph2, bx2, pgn5, pgn10, pgn100, rp8, rp60, )
# make 3D shapes
sph3 = make_sphere3(rand(5))
bx3 = make_box3(rand(16))
ppr6 = make_6polyprism(rand(16))
cylpr = make_cyl(rand(9))
shapes3D = ( sph3, bx3, ppr6, cylpr )

@testset verbose = true "surfpt_nearby AD Gradients" begin
    @testset "2D x-vector gradients of surfpt_nearby(x,shape) for 2D shape type $(typeof(sh))" for sh in shapes2D
        p = rand(2)
        @test test_AD(x->sum(sum(surfpt_nearby(SVector{2}(x),sh))),p)
    end
    @testset "3D x-vector gradients of surfpt_nearby(x,shape) for 2D shape type $(typeof(sh))" for sh in shapes2D
        p = rand(3)
        @test test_AD(x->sum(sum(surfpt_nearby(SVector{3}(x),sh))),p)
    end
    @testset "3D x-vector gradients of surfpt_nearby(x,shape) for 3D shape type $(typeof(sh))" for sh in shapes3D
        p = rand(3)
        @test test_AD(x->sum(sum(surfpt_nearby(SVector{3}(x),sh))),p)
    end
    @testset "2D shape parameter gradients of surfpt_nearby(x::SVector{2},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
        xy = SVector{2}(rand(2))
        p = rand(np)
        @test test_AD(x->sum(sum(surfpt_nearby(xy,fsh(x)))),p)
    end
    @testset "2D shape parameter gradients of surfpt_nearby(x::SVector{3},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh2D,np_fsh2D)
        xyz = SVector{3}(rand(3))
        p = rand(np)
        @test test_AD(x->sum(sum(surfpt_nearby(xyz,fsh(x)))),p)
    end
    @testset "3D shape parameter gradients of surfpt_nearby(x::SVector{3},shape) for shape fn $fsh with $np params" for (fsh,np) in zip(fsh3D,np_fsh3D)
        xyz = SVector{3}(rand(3))
        p = rand(np)
        @test test_AD(x->sum(sum(surfpt_nearby(xyz,fsh(x)))),p)
    end
end


@test test_AD(f_spn_22,p0)

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

##
Zygote.gradient()

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
