@testset "triangular cylinder 3D" begin
    vxl = (SVector(0,0,0), SVector(1,1,1))
    nout = SVector(1,1,0)

    @test_nowarn (r₀ = SVector(0.5,0,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVector(0.5,0,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVector(0.5,0,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 3D"

@testset "triangular cylinder 2D" begin
    vxl = (SVector(0,0), SVector(1,1))
    nout = SVector(1,1)

    @test_nowarn (r₀ = SVector(0.5,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVector(0.5,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVector(0.5,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVector(1,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 2D"

@testset "quadrangular cylinder 3D" begin
    @test_nowarn @inferred(volfrac((SVector(0,0,0),SVector(1,1,1)), SVector(1,2,0), SVector(0.5,0.5,0.5)))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = mean(vxl)
            nout = randn(3)
            nout[rand(1:3)] = 0
            result &= volfrac(vxl, SVector{3}(nout), r₀)≈0.5
        end
        result
    end
end  # @testset "quadrangular cylinder 3D"

@testset "quadrangular cylinder 2D" begin
    @test_nowarn @inferred(volfrac((SVector(0,0),SVector(1,1)), SVector(1,2), SVector(0.5,0.5)))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(2)), @SVector(rand(2)))
            r₀ = mean(vxl)
            nout = @SVector randn(2)
            result &= volfrac(vxl, nout, r₀)≈0.5
        end
        result
    end
end  # @testset "quadrangular cylinder 2D"

@testset "general cases" begin
    # Test random cases.
    @test_nowarn @inferred(volfrac((-@SVector(rand(3)),@SVector(rand(3))), SVector(0,0,0), @SVector(randn(3))))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = mean(vxl)
            nout = @SVector randn(3)
            result &= volfrac(vxl, nout, r₀)≈0.5
        end
        result
    end

    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = @SVector randn(3)
            nout = @SVector randn(3)
            result &= (volfrac(vxl, nout, r₀) + volfrac(vxl, -nout, r₀) ≈ 1)
        end
        result
    end

    # Test boundary cases.
    vxl = (SVector(0,0,0), SVector(1,1,1))
    r₀ = SVector(0,0,0)
    @test (nout = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,0,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(1,1,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,-1,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(1,1,1); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,-1,-1); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(-1,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_tricyl()
    @test (nout = SVector(-1,-1,1); volfrac(vxl, nout, r₀) ≈ 5/6)  # rvol_gensect()
    @test (nout = SVector(1,-2,1); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadsect()
    r₀ = SVector(0.5,0.5,0)
    @test (nout = SVector(-2,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

    # Test tolerance to floating-point arithmetic.
    vxl = (SVector(0,9.5,0), SVector(1,10.5,1))
    r₀ = SVector(0.5,10.0,0.5)
    nout = SVector(0.0,1/√2,1/√2)
    @test volfrac(vxl, nout, r₀) ≈ 0.5

    # Test rvol_quadsect() for nontrivial cases.
    vxl = (SVector(0,0,0), SVector(1,1,2))
    r₀ = SVector(0.5, 0.5, 0.5)
    @test begin
        result = true
        for i = 1:100
            nout = SVector(randn()/20, randn()/20, 1)
            result &= volfrac(vxl, nout, r₀)≈0.5/2
        end
        result
    end
end  # @testset "general cases"

@testset "non-SVector arguments" begin
    # `volfrac` must accept plain `AbstractVector`s (e.g. an `MVector`), not just `SVector`s.
    # Some AD tools (notably Zygote) hand `volfrac` the outward normal returned by
    # `surfpt_nearby` as an `MVector` during the forward pass, which previously hit a
    # MethodError; the fallback method converts the arguments to `SVector`s.
    vxl3 = (SVector(-0.6,-0.7,-0.3), SVector(1.4,1.3,1.7))
    nout3 = SVector(0.5,0.6,0.7)
    r₀3 = SVector(0.5,0.5,0.9)
    ref3 = volfrac(vxl3, nout3, r₀3)

    # Non-SVector normal (the reported case): MVector and Vector must match the SVector call.
    @test volfrac(vxl3, MVector(0.5,0.6,0.7), r₀3) == ref3
    @test volfrac(vxl3, [0.5,0.6,0.7], r₀3) == ref3
    # Every argument as a non-SVector vector.
    @test volfrac(([-0.6,-0.7,-0.3], [1.4,1.3,1.7]), [0.5,0.6,0.7], [0.5,0.5,0.9]) == ref3
    # Mixed: SVector voxel, MVector normal, Vector r₀.
    @test volfrac(vxl3, MVector(0.5,0.6,0.7), [0.5,0.5,0.9]) == ref3

    # 2D and 1D likewise.
    vxl2 = (SVector(-0.6,-0.7), SVector(1.4,1.3))
    ref2 = volfrac(vxl2, SVector(0.5,0.6), SVector(0.5,0.5))
    @test volfrac(vxl2, MVector(0.5,0.6), SVector(0.5,0.5)) == ref2
    vxl1 = (SVector(-0.6,), SVector(1.4,))
    ref1 = volfrac(vxl1, SVector(0.7,), SVector(0.2,))
    @test volfrac(vxl1, [0.7], [0.2]) == ref1

    # Mismatched lengths are rejected with a clear error.
    @test_throws DimensionMismatch volfrac(vxl3, [0.5,0.6], r₀3)
end  # @testset "non-SVector arguments"
