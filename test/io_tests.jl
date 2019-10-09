using NeRCA
using Test


const DETX = joinpath(@__DIR__, "data", "detx_v3.detx")
const H5FILE = joinpath(@__DIR__, "data", "mupage.h5")


@testset "calibration" begin
    calib = NeRCA.read_calibration(DETX)

    @test 23 == calib.det_id
    @test 6 == length(values(calib.pos))
    @test 6 == length(values(calib.dir))
    @test 6 == length(values(calib.t0))
    @test 6 == length(values(calib.du))
    @test 6 == length(values(calib.floor))
    @test 18 == length(values(calib.omkeys))
    @test 23.9 ≈ calib.max_z
    @test 2 == calib.n_dus
    @test 3.6 ≈ calib.pos[3][2].z
end


@testset "reading KM3HDF5 hits" begin
    hits = read_hits(H5FILE, 0)
    @test 4886 == length(hits)
    @test 3.3263139e7 ≈ hits[end].t
    hits = read_hits(H5FILE, 2)
    @test 4219 == length(hits)
    @test 5.1047085e7 ≈ hits[end].t
end


@testset "reading KM3HDF5 MC hits" begin
    hits = NeRCA.read_mchits(H5FILE, 0)
    @test 399 == length(hits)
    @test 4429.716844306684 ≈ hits[end].t
end


@testset "reading KM3HDF5 MC tracks" begin
    mc_tracks = NeRCA.read_mctracks(H5FILE)
    @test 3 == length(mc_tracks)
    @test 11 == length(mc_tracks[0])
    @test 3 == length(mc_tracks[2])
end
