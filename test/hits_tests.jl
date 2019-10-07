using NeRCA
using Test


hits = [Hit(1, 8, 100, 20, false),
        Hit(2, 9, 101, 21, true),
        Hit(3, 8, 112, 22, true),
        Hit(4, 8, 114, 23, false),
        Hit(5, 8, 134, 24, true),
        Hit(6, 8, 156, 25, false),
        Hit(7, 8, 133, 26, true),
        Hit(8, 8, 145, 26, false),
        Hit(9, 10, 900, 26, false),
        Hit(1, 11, 1000, 26, false),
        Hit(1, 11, 1001, 26, false),
        Hit(1, 11, 1002, 26, false)]


# triggered()
thits = triggered(hits)
@test 4 == length(thits)
@test 9 == thits[1].dom_id


# nfoldhits()
twofoldhits = nfoldhits(hits, 10, 2)
@test 4 == length(twofoldhits)
threefoldhits = nfoldhits(hits, 15, 3)
@test 3 == length(threefoldhits)


# domhits()
dhits = domhits(hits)
@test 7 == length(dhits[8])
@test 20 == dhits[8][1].tot
@test dhits[8][6].triggered

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
