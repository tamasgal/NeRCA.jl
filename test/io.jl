using KM3NeT
using Base.Test

hits = [RawHit(1, 8, 100, 20, false),
        RawHit(2, 9, 101, 21, true),
        RawHit(3, 8, 112, 22, true),
        RawHit(4, 8, 114, 23, false),
        RawHit(5, 8, 134, 24, true),
        RawHit(6, 8, 156, 25, false),
        RawHit(7, 8, 133, 26, true),
        RawHit(8, 8, 145, 26, false)]

println("triggered()")
thits = triggered(hits)
@test 4 == length(thits)
@test 9 == thits[1].dom_id

println("nfoldhits()")
twofoldhits = nfoldhits(hits, 10, 2)
@test 4 == length(twofoldhits)
threefoldhits = nfoldhits(hits, 15, 3)
@test 3 == length(threefoldhits)