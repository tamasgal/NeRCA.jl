using NeRCA
using KM3io
using Test

const ONLINEFILE = joinpath(@__DIR__, "data", "km3net_online.root")

hits = [TriggeredHit(8, 1, 100, 20, false),
        TriggeredHit(9, 2, 101, 21, true),
        TriggeredHit(8, 3, 112, 22, true),
        TriggeredHit(8, 4, 114, 23, false),
        TriggeredHit(8, 5, 134, 24, true),
        TriggeredHit(8, 6, 156, 25, false),
        TriggeredHit(8, 7, 133, 26, true),
        TriggeredHit(8, 8, 145, 26, false),
        TriggeredHit(10, 9, 900, 26, false),
        TriggeredHit(11, 1, 1000, 26, false),
        TriggeredHit(11, 1, 1001, 26, false),
        TriggeredHit(11, 1, 1002, 26, false)]


# triggered()
@testset "triggered()" begin
    thits = triggered(hits)
    @test 4 == length(thits)
    @test 9 == thits[1].dom_id
end

# nfoldhits()
@testset "nfoldhits()" begin
    twofoldhits = nfoldhits(hits, 10, 2)
    @test 4 == length(twofoldhits)
    threefoldhits = nfoldhits(hits, 15, 3)
    @test 3 == length(threefoldhits)
end

# domhits()
@testset "domhits()" begin
    dhits = domhits(hits)
    @test 7 == length(dhits[8])
    @test 20 == dhits[8][1].tot
    @test 1 == dhits[8][6].trigger_mask
end

# multiplicities
sorted_hits = sort(hits, by=h->h.t)
sort!(sorted_hits, by=h->h.dom_id)
mtps, mtp_ids = NeRCA.count_multiplicities(sorted_hits, 10)
@test (1, 1) == (mtps[1], mtp_ids[1])
@test (2, 2) == (mtps[2], mtp_ids[2])
@test (2, 2) == (mtps[3], mtp_ids[3])
@test (2, 3) == (mtps[4], mtp_ids[4])
@test (2, 3) == (mtps[5], mtp_ids[5])
@test (1, 4) == (mtps[6], mtp_ids[6])
@test (1, 5) == (mtps[7], mtp_ids[7])
@test (1, 6) == (mtps[8], mtp_ids[8])
@test (1, 7) == (mtps[9], mtp_ids[9])
@test (3, 8) == (mtps[10], mtp_ids[10])
@test (3, 8) == (mtps[11], mtp_ids[11])
@test (3, 8) == (mtps[12], mtp_ids[12])

# combine
f = OnlineFile(ONLINEFILE)
for event in f.events
    local hits = NeRCA.combine(event.snapshot_hits, event.triggered_hits)
    @test length(hits) == length(event.snapshot_hits)
    @test length(NeRCA.triggered(hits)) == length(event.triggered_hits)
end
