# Even though shapes are specified with real geometric parameters, the size of voxels can be
# complex when the underlying space is stretched by complex factors.  So, we allow complex
# `vxl` below.  Even in such a case, r₀ and nout are still real.  (Maybe this condition
# should be relaxed later.)

export volfrac

using Zygote: ignore, @ignore, dropgrad

const X, Y, Z = 1, 2, 3
const XYZ, YZX, ZXY = (X,Y,Z), (Y,Z,X), (Z,X,Y)
const UVW = YZX, ZXY, XYZ
const N, P = 1, 2  # negative, positive
const NP = (N, P)

corner(vxl::NTuple{2,SVector{3,<:Number}}, sx::Integer, sy::Integer, sz::Integer) =
    SVector(vxl[sx][X], vxl[sy][Y], vxl[sz][Z])

# Calculate the bit array that indicates corner contained-ness.
#
# The bit array is 8 bits (= 1 byte).  The kth bit is 1 if the kth corner is contained in
# the half space defined by the plane, an 0 otherwise.  The 8 corners of the voxel are
# indexed 1 through 8 in the order of (---), (+--), (-+-), (++-), (--+), (+-+), (-++), (+++),
# where, e.g., (+--) indicates the corner at the +x, -y, -z corner.
function corner_bits(vxl::NTuple{2,SVector{3,<:Number}},  # two ends of solid diagonal of voxel
                     nout::SVector{3,<:Real}, # unit outward normal of plane
                     nr₀::Real)  # equation of plane: nout⋅(r - r₀) = 0, or nout⋅r = nr₀
    cbits = 0x00
    bit = 0x01
    n_on = 0
    for sz = NP, sy = NP, sx = NP
        r = corner(vxl, sx, sy, sz)
        nr = nout⋅r
        if nr ≤ nr₀  # corner is contained and not on boundary (nout is outward normal)
            cbits |= bit
            n_on += nr==nr₀  # corner is on boundary
        end
        bit <<= 1
    end

    return cbits, n_on
end

# Determine if the corner contained-ness bit array corresponds to the case where the plane
# crosses one of the three sets of four parallel edges of the voxel.
function isquadsect(cbits::UInt8)
    return cbits==0x0F || cbits==~0x0F || cbits==0x33 || cbits==~0x33 || cbits==0x55 || cbits==~0x55

    # Equivalent to
    # n = count_ones(cbits)
    # m = count_ones(cbits & 0x0F)
    # return n==4 && iseven(m)
end

# For the cases where the plane crosses a set of four parallel edges of the voxel, determine
# which direction those edges lie.
function edgedir_quadsect(cbits::UInt8)
    if cbits==0x0F || cbits==~0x0F
        dir = Z
    elseif cbits==0x33 || cbits==~0x33
        dir = Y
    else
        @assert cbits==0x55 || cbits==~0x55
        dir = X
    end

    return dir
end

# Return the volume fraction when the plane croses a set of four parallel edges.
function rvol_quadsect(vxl::NTuple{2,SVector{3,<:Number}}, nout::SVector{3,<:Real}, nr₀, cbits::UInt8)
    w = edgedir_quadsect(cbits)
    ∆w = vxl[P][w] - vxl[N][w]

    u, v, _w = UVW[w]
    nu, nv, nw = nout[u], nout[v], nout[w]
    mean_cepts = 4nr₀
    for sv = NP, su = NP
        mean_cepts -= nu*vxl[su][u] + nv*vxl[sv][v]
    end
    mean_cepts /=  nw * 4∆w

    sw = nw>0 ? N : P  # nw cannot be 0
    return abs(mean_cepts - vxl[sw][w]/∆w)
end

# Calculate the volume fraction for the most general cases, by cutting out corners from a
# triangular pyramid.
#
# Assume count_ones(cbits) ≤ 4.  Othewise, call this function with flipped nout, nr₀,
# cbits.
function rvol_gensect(vxl::NTuple{2,SVector{3,<:Number}}, nout::SVector{3,<:Real}, nr₀::Real, cbits::UInt8)
    s = (nout.<0) .+ 1
    c = corner(vxl, s[X], s[Y], s[Z])  # corner coordinates
    ∆ = vxl[P] - vxl[N]  # vxl edges
    nc = nout .* c
    #rmax0, rmid0, rmin0 =  abs.(((nr₀-sum(nc)) .+ nc) ./ nout - c) ./ ∆ # (lengths from corner to intercetps) / (voxel edges)
    nc_sum = nc[1] + nc[2] + nc[3]
    # rs =  (((nr0-nc_sum) .+ nc) ./ nout - c) ./ ∆
    rmax0 =  abs(((nr₀-nc_sum) + nc[1]) / nout[1] - c[1]) / ∆[1]
    rmid0 =  abs(((nr₀-nc_sum) + nc[2]) / nout[2] - c[2]) / ∆[2]
    rmin0 =  abs(((nr₀-nc_sum) + nc[3]) / nout[3] - c[3]) / ∆[3]

    # println("####### inside rvol_gensect ########")
    # println("s: $(s), [$(typeof(s))]")
    # println("c: $(c), [$(typeof(c))]")
    # println("Δ: $(∆), [$(typeof(∆))]")
    # println("nc: $(nc), [$(typeof(nc))]")
    # println("### rvals, presort ####")
    # println("rmax0: $(rmax0), [$(typeof(rmax0))]")
    # println("rmid0: $(rmid0), [$(typeof(rmid0))]")
    # println("rmin0: $(rmin0), [$(typeof(rmin0))]")

    # nx, ny, nz = nout
    # sx, sy, sz = ((nx≥0 ? N : P), (ny≥0 ? N : P), (nz≥0 ? N : P))  # signs of corner
    # cx, cy, cz = vxl[nX][sx], vxl[nY][sy], vxl[nZ][sz]  # corner coordinates
    # ∆x, ∆y, ∆z = vxl[nX][nP]-vxl[nX][nN], vxl[nY][nP]-vxl[nY][nN], vxl[nZ][nP]-vxl[nZ][nN]  # voxel edges
    # nxcx, nycy, nzcz = nx*cx, ny*cy, nz*cz
    # rmax, rmid, rmin =  # (lengths from corner to intercetps) / (voxel edges)
    #     abs((nr₀-nycy-nzcz)/nx-cx)/∆x, abs((nr₀-nzcz-nxcx)/ny-cy)/∆y, abs((nr₀-nxcx-nycy)/nz-cz)/∆z

    # Sort rmax, rmin, rmin properly.
    # if rmax < rmid; rmax, rmid = rmid, rmax; end
    # if rmid < rmin; rmid, rmin = rmin, rmid; end
    # if rmax < rmid; rmax, rmid = rmid, rmax; end

    if rmax0 < rmid0
        rmax1 = rmid0
        rmid1 = rmax0
        rmin1 = rmin0
    else
        rmax1 = rmax0
        rmid1 = rmid0
        rmin1 = rmin0
    end

    if rmid1 < rmin1
        rmid2 = rmin1
        rmin2 = rmid1
        rmax2 = rmax1
    else
        rmid2 = rmid1
        rmin2 = rmin1
        rmax2 = rmax1
    end

    if rmax2 < rmid2
        rmax = rmid2
        rmid = rmax2
        rmin = rmin2
    else
        rmax = rmax2
        rmid = rmid2
        rmin = rmin2
    end

    # println("### rvals, post-sort ####")
    # println("rmax: $(rmax), [$(typeof(rmax))]")
    # println("rmid: $(rmid), [$(typeof(rmid))]")
    # println("rmin: $(rmin), [$(typeof(rmin))]")

    # Calculate the volume of the triangular pyramid, and cut off appropriate corners.
    # tmax = 1 - 1/rmax
    # rvol_core = 1 + tmax + tmax^2

    # println("### tmax,rvol_core vals before check ####")
    # println("tmax: $(tmax), [$(typeof(tmax))]")
    # println("rvol_core: $(rvol_core), [$(typeof(rvol_core))]")

    # Below, if rmax == Inf, rmid and rmin must be ≤ 1 in exact arithmetic and subtraction
    # must not occur, but they can be > 1 in floating-point arithmetic and lead to large
    # subtraction.  Prevent this by subtracting only for rmax ≠ Inf.
    # isfin_rmax = !isinf(rmax)
    # if rmid > 1 && isfin_rmax
    #     tmid = 1 - 1/rmid
    #     rvol_core -= rmax * tmid^3
    # end
    # if rmin > 1 && isfin_rmax
    #     tmin = 1 - 1/rmin
    #     rvol_core -= rmax * tmin^3
    # end


    if !isinf(rmax)
        tmax = 1 - 1/rmax
        rvol_core0 = 1 + tmax + tmax^2
        if rmid > 1
            tmid = 1 - 1/rmid
            rvol_core1 = rvol_core0 - ( rmax * tmid^3 )
        else
            rvol_core1 = rvol_core0
        end

        if rmin > 1
            tmin = 1 - 1/rmin
            rvol_core = rvol_core1 - ( rmax * tmin^3 )
        else
            rvol_core = rvol_core1
        end
    else
        tmax = 1. # 1 - 1/rmax
        rvol_core = 3. # 1 + tmax + tmax^2
    end

    # println("### tmax,rvol_core vals after check ####")
    # println("tmax: $(tmax), [$(typeof(tmax))]")
    # println("rvol_core: $(rvol_core), [$(typeof(rvol_core))]")

    return rvol_core * rmin * rmid / 6.
    # out = rvol_core * rmin * rmid / 6.
    # out = @Zygote.showgrad(rvol_core) * @Zygote.showgrad(rmin) * @Zygote.showgrad(rmid) / 6.
    # return @Zygote.showgrad(out)
end

"""
    volfrac(voxel, nout, r₀)

Returns the volume fraction `rvol = vol(voxel ⋂ half-space) / vol(voxel)` that indicates how
much portion of the volume of a voxel is included the half-space defined by a plane.  The
result is between 0.0 and 1.0, inclusive.  Inside and outside of the half space is
determined by `nout`, which is the outward normal.  (If a point is on the `nout` side of the
plane, then it is outside.)

The voxel needs to be aligned with the Cartesian axes, and is described by `v =
([xn,yn,zn], [xp,yp,zp])` that specifies the x-, y-, z-boundaries of the voxel.  The
half-space is described by the boundary plane.  The boundary plane is described by the
outward normal vector `nout = [nx, ny, nz]` and a point `r₀ = [rx, ry, rz]` on the plane.
`nout` does not have to be normalized.
"""
function volfrac(vxl::NTuple{2,SVector{3,<:Number}}, nout::SVector{3,<:Real}, r₀::SVector{3,<:Real})
    nr₀ = nout⋅r₀
    cbits, n_on = corner_bits(vxl, nout, nr₀)
    n_in = count_ones(cbits)  # number of corners contained

    # nr₀ = dot(@Zygote.showgrad(nout),@Zygote.showgrad(r₀))
    # cbits, n_on = corner_bits(@Zygote.showgrad(vxl), @Zygote.showgrad(nout), @Zygote.showgrad(nr₀))

    # nr₀ = dot(@show(nout),@show(r₀))
    # cbits, n_on = corner_bits(@show(vxl), @show(nout), @show(nr₀))
    # n_in = count_ones(@show(cbits))  # number of corners contained
    # println("vxl: $(vxl), [$(typeof(vxl))]")
    # println("nout: $(nout), [$(typeof(nout))]")
    # println("r₀: $(r₀), [$(typeof(r₀))]")

    # println("nr₀: $(nr₀), [$(typeof(nr₀))]")
    # println("cbits: $(cbits), [$(typeof(cbits))]")
    # println("n_on: $(n_on), [$(typeof(n_on))]")
    # println("n_in: $(n_in), [$(typeof(n_in))]")


    if n_in == 8  # voxel is inside half-space
        # println("path 1: voxel is inside half-space")
        rvol = 1.
    elseif n_in - n_on == 0  # voxel is outside half-space
        # println("path 2: voxel is outside half-space")
        rvol = 0.
    elseif isquadsect(cbits) # plane crosses a set of four parallel edges of voxel
        # println("path 3: plane crosses a set of four parallel edges of voxel")
        rvol = rvol_quadsect(vxl, nout, nr₀, cbits)
        # rvol = rvol_quadsect(@Zygote.showgrad(vxl), @Zygote.showgrad(nout), @Zygote.showgrad(nr₀), @Zygote.showgrad(cbits))
    elseif n_in ≤ 4 # general cases with n_in ≤ 4
        # println("path 4: general cases with n_in ≤ 4")
        rvol = rvol_gensect(vxl, nout, nr₀, cbits)
        # rvol = rvol_gensect(@Zygote.showgrad(vxl), @Zygote.showgrad(nout), @Zygote.showgrad(nr₀), @Zygote.showgrad(cbits))
    else  # general cases with n_in ≥ 5
        # println("path 5: general cases with n_in ≥ 5")
        @assert n_in ≥ 5
        rvol = 1. - rvol_gensect(vxl, .-nout, -nr₀, ~cbits)
    end

    return rvol
end

volfrac(vxl::NTuple{2,SVector{2,<:Number}}, nout::SVector{2,<:Real}, r₀::SVector{2,<:Real}) =
    volfrac((SVector(vxl[N][1],vxl[N][2],0), SVector(vxl[P][1],vxl[P][2],1)),
            SVector(nout[1], nout[2], 0),
            SVector(r₀[1], r₀[2], 0))

volfrac(vxl::NTuple{2,AbstractVector{<:Real}}, nout::AbstractVector{<:Real}, r₀::AbstractVector{<:Real}) =
    volfrac((SVector(vxl[N][1],vxl[N][2],0), SVector(vxl[P][1],vxl[P][2],1)),
            SVector(nout[1], nout[2], 0),
            SVector(r₀[1], r₀[2], 0))
