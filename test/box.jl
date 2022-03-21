@testset "Box" begin
    b = Box([0,0], [2,4])  # specify center and radii
    @test b == deepcopy(b)
    @test hash(b) == hash(deepcopy(b))
    @test [0.3,-1.5] ∈ b
    @test [0.3,-2.5] ∉ b

    @test ((x,nout) = surfpt_nearby([0,0],b); (x≈[1,0] && nout≈[1,0]) || (x≈[0,2] && nout≈[0,1]))  # handle point at center properly
    @test all([(p = [1sx,2sy]; surfpt_nearby(1.1p,b) == (p, normalize(1.1p-p))) for sx = (-1,1), sy = (-1,1)])  # outside corners
    @test all([((x, nout) = surfpt_nearby([1sx,2sy],b); x≈[1sx,2sy] && all([sx 0; 0 sy]*nout.≥0) && norm(nout)≈1) for sx = (-1,1), sy = (-1,1)])  # at corners
    @test all([surfpt_nearby([ρ*1sx,1sy],b) == ([1sx,1sy], [sx,0]) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,1), sy = (-1,0,1)])  # around faces
    @test all([surfpt_nearby([0.5sx,ρ*2sy],b) == ([0.5sx,2sy], [0,sy]) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,0,1), sy = (-1,1)])  # around faces

    @test normal([1.1,0],b) == [1,0]
    @test normal([-1.1,0],b) == [-1,0]
    @test normal([1.1,2.1],b) == [1,1]/√2
    @test bounds(b) == ([-1,-2],[1,2])
    # @test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
    @test bounds(Box(SVector{2,Float64}([0,0]), SVector{2,Float64}([2,4]), SMatrix{2,2,Float64}([1 1; 1 -1]))) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
    @test checkbounds(b)
    # @test checkbounds(Box([0,0], [2,4], [1 1; 1 -1]))
    @test checkbounds(Box(SVector{2,Float64}([0,0]), SVector{2,Float64}([2,4]), SMatrix{2,2,Float64}([1 1; 1 -1])))
    # @test (∆ = rand(2); translate(b,∆) ≈ Box([0,0]+∆, [2,4]))
    @test (∆ = rand(2); translate(b,∆) ≈ Box(SVector{2,Float64}([0,0]+∆), SVector{2,Float64}([2,4])) )
end  # @testset "Box"


@testset "Box, rotated" begin
    # ax1, ax2 = normalize.(([1,-1], [1,1]))
    # r1, r2 = 1, 2  # "radii"
    # br = Box([0,0], [2r1, 2r2], [ax1 ax2])
    # R = [ax1 ax2]  # rotation matrix
    ax1, ax2 = normalize.((SVector{2,Float64}([1,-1]), SVector{2,Float64}([1,1])))
    r1, r2 = 1.0, 2.0  # "radii"
    br = Box(SVector{2,Float64}([0,0]), SVector{2,Float64}([2r1, 2r2]), SMatrix{2,2,Float64}([ax1 ax2]))
    R = SMatrix{2,2,Float64}([ax1 ax2])
    Cin = R * (GeometryPrimitives.signmatrix(br) .* (one⁻ .* [r1,r2]))  # around corners, inside
    Cout = R * (GeometryPrimitives.signmatrix(br) .* (one⁺ .* [r1,r2]))  # around corners, outside
    for j = 1:2; @test Cin[:,j] ∈ br; end
    for j = 1:2; @test Cout[:,j] ∉ br; end

    @test br == deepcopy(br)
    @test hash(br) == hash(deepcopy(br))

    @test ((x,nout) = surfpt_nearby([0,0],br); (x≈r1*ax1 && nout≈ax1) || (x≈r2*ax2 && nout≈ax2))  # handle point at center properly
    @test all([(p = (s1*r1*ax1+s2*r2*ax2); surfpt_nearby(1.1p,br) ≈ (p,normalize(1.1p-p))) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
    @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,br); x≈p && all([s1*ax1 s2*ax2]'*nout.≥0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
    @test all([(p1 = s1*r1*ax1; p2 = s2*r2/2*ax2; surfpt_nearby(ρ*p1+p2,br) ≈ (p1+p2, s1*ax1)) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,1), s2 = (-1,0,1)])  # around faces
    @test all([(p1 = s1*r1/2*ax1; p2 = s2*r2*ax2; surfpt_nearby(p1+ρ*p2,br) ≈ (p1+p2, s2*ax2)) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,0,1), s2 = (-1,1)])  # around faces

    @test normal(R*[1.1r1, 0], br) ≈ R*[1,0]
    @test normal(R*[-1.1r1, 0], br) ≈ R*[-1,0]
    @test normal(R*[0, 1.1r2], br) ≈ R*[0,1]
    @test normal(R*[0, -1.1r2], br) ≈ R*[0,-1]
    @test normal(R*[1.1r1, 1.01r2], br) ≈ R*(normalize([1.1r1, 1.01r2]-[r1,r2]))

    xmax = (R*[r1,r2])[1]
    ymax = (R*[-r1,r2])[2]
    @test bounds(br) ≈ (-[xmax,ymax], [xmax,ymax])
    @test checkbounds(br)

    @test (∆ = rand(2); translate(br,∆) ≈ Box([0,0]+∆, [2r1, 2r2], [ax1 ax2]))
end  # @testset "Box, rotated"


@testset "Box, skewed" begin
    ax1, ax2 = normalize.(([1,-1], [0,1]))
    r1, r2 = 1, 1  # "radii"
    bs = Box([0,0], [2r1, 2r2], [ax1 ax2])

    @test bs == deepcopy(bs)
    @test hash(bs) == hash(deepcopy(bs))

    n1, n2 = normalize.(([1,0], [1,1]))
    @test ((x,nout) = surfpt_nearby([0,0],bs); (x⋅n1≈1/√2 && nout≈n1) || (x⋅n2≈1/√2 && nout≈n2))  # handle point at center properly
    @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(1.1p,bs); nout⋅(s1*n1+s2*n2)>0 && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
    @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,bs); x≈p && nout⋅(s1*n1+s2*n2)>0 && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
    # Use the above, less demanding tests instead of the below two; surfpt_nearby suffers from some inaccuracy around corners of skewed boxes.
    # @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(1.1p,bs); all([s1*n1 s2*n2]'*nout.>0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
    # @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,bs); all([s1*n1 s2*n2]'*nout.≥0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
    @test all([(p1 = s1*r1*ax1; p2 = s2*r2/2*ax2; surfpt_nearby(ρ*p1+p2,bs) ≈ (p1+p2, s1*n1)) for ρ = (one⁻,1,one⁺), s1 = (-1,1), s2 = (-1,0,1)])  # around faces
    @test all([(p1 = s1*r1/2*ax1; p2 = s2*r2*ax2; surfpt_nearby(p1+ρ*p2,bs) ≈ (p1+p2, s2*n2)) for ρ = (one⁻,1,one⁺), s1 = (-1,0,1), s2 = (-1,1)])  # around faces

    @test norm(normal([0,1], bs)) ≈ 1

    xmax = (r1*ax1+r2*ax2)[1]
    ymax = (r2*ax2-r1*ax1)[2]
    @test bounds(bs) ≈ (-[xmax,ymax],[xmax,ymax])
    @test checkbounds(bs)

    @test (∆ = rand(2); translate(bs,∆) ≈ Box([0,0]+∆, [2r1, 2r2], [ax1 ax2]))
end  # @testset "Box, skewed"
