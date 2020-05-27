using NeRCA
using Test


const DETX = joinpath(@__DIR__, "data", "detx_v3.detx")
const H5FILE = joinpath(@__DIR__, "data", "mupage.h5")
const ONLINEFILE = joinpath(@__DIR__, "data", "km3net_online.root")
const OFFLINEFILE = joinpath(@__DIR__, "data", "km3net_offline.root")


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

@testset "km3net online files" begin
    f = NeRCA.OnlineFile(ONLINEFILE)
    hits = NeRCA.read_hits(f)
    @test 3 == length(hits)  # grouped by event
    @test 96 == length(hits[1])
    @test [806451572, 806451572, 806455814] == [h.dom_id for h in hits[1][1:3]]
    @test [809524432, 809524432, 809544061] == [h.dom_id for h in hits[1][end-2:end]]
    @test [30733918, 30733916, 30733256] == [h.time for h in hits[1][1:3]]
    @test [30733864, 30734686, 30735112] == [h.time for h in hits[1][end-2:end]]
    @test [10, 13, 0, 3, 1] == [h.channel_id for h in hits[1][1:5]]
    @test [10, 10, 22, 24, 17] == [h.channel_id for h in hits[1][end-4:end]]
    @test [26, 19, 25, 22, 28] == [h.tot for h in hits[1][1:5]]
    @test [6, 10, 29, 28, 27] == [h.tot for h in hits[1][end-4:end]]
    @test 124 == length(hits[2])
    @test [806455814, 806483369, 806483369] == [h.dom_id for h in hits[2][1:3]]
    @test [809521500, 809526097, 809526097, 809544058, 809544061] == [h.dom_id for h in hits[2][end-4:end]]
    @test [58728018, 58728107, 58729094] == [h.time for h in hits[2][1:3]]
    @test [58729410, 58729741, 58729262] == [h.time for h in hits[2][end-4:end]]
    @test [15, 5, 14, 23, 9] == [h.channel_id for h in hits[2][1:5]]
    @test [17,  5, 18, 24,  8] == [h.channel_id for h in hits[2][end-4:end]]
    @test [27, 24, 21, 17, 22] == [h.tot for h in hits[2][1:5]]
    @test [21, 23, 25, 27, 27] == [h.tot for h in hits[2][end-4:end]]
    @test 78 == length(hits[3])
    @test [806451572, 806483369, 806483369] == [h.dom_id for h in hits[3][1:3]]
    @test [809526097, 809526097, 809526097, 809544058, 809544061] == [h.dom_id for h in hits[3][end-4:end]]
    @test [63512204, 63511134, 63512493] == [h.time for h in hits[3][1:3]]
    @test [809526097, 809526097, 809526097] == [h.time for h in hits[3][end-4:end]]
    @test [4, 9, 5, 17, 20] == [h.channel_id for h in hits[3][1:5]]
    @test [5,  7, 24, 23, 10] == [h.channel_id for h in hits[3][end-4:end]]
    @test [26, 29, 30, 23, 30] == [h.tot for h in hits[3][1:5]]
    @test [28, 11, 27, 24, 23] == [h.tot for h in hits[3][end-4:end]]
end


# @testset "reading KM3HDF5 hits" begin
#     hits = read_hits(H5FILE, 0)
#     @test 4886 == length(hits)
#     @test 3.3263139e7 ≈ hits[end].t
#     hits = read_hits(H5FILE, 2)
#     @test 4219 == length(hits)
#     @test 5.1047085e7 ≈ hits[end].t
# end


# @testset "reading KM3HDF5 MC hits" begin
#     hits = NeRCA.read_mchits(H5FILE, 0)
#     @test 399 == length(hits)
#     @test 4429.716844306684 ≈ hits[end].t
# end


# @testset "reading KM3HDF5 MC tracks" begin
#     mc_tracks = NeRCA.read_mctracks(H5FILE)
#     @test 3 == length(mc_tracks)
#     @test 11 == length(mc_tracks[0])
#     @test 3 == length(mc_tracks[2])
# end
