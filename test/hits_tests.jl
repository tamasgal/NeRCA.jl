using NeRCA
using KM3io
using KM3NeTTestData
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


# nfoldhits()
@testset "nfoldhits()" begin
    twofoldhits = nfoldhits(hits, 10, 2)
    @test 4 == length(twofoldhits)
    threefoldhits = nfoldhits(hits, 15, 3)
    @test 3 == length(threefoldhits)
end

@testset "modulemap()" begin
    dhits = modulemap(hits)
    @test 7 == length(dhits[8])
    @test 20 == dhits[8][1].tot
    @test 1 == dhits[8][6].trigger_mask
end

@testset "sort hits" begin
    rhits = [
        HitR1(1, zero(Position), 15, 0, 1, 1),
        HitR1(2, zero(Position), 14, 0, 1, 1),
        HitR1(2, zero(Position), 11, 0, 1, 1),
        HitR1(2, zero(Position), 10, 0, 1, 1),
        HitR1(3, zero(Position), 12, 0, 1, 1),
        HitR1(2, zero(Position), 14, 0, 1, 1),
        HitR1(1, zero(Position), 13, 0, 1, 1),
    ]

    sort!(rhits)

    @test 1 == rhits[1].dom_id
    @test 1 == rhits[2].dom_id
    @test 2 == rhits[3].dom_id
    @test 2 == rhits[4].dom_id
    @test 2 == rhits[5].dom_id
    @test 2 == rhits[6].dom_id
    @test 3 == rhits[7].dom_id

    @test 13 == rhits[1].t
    @test 15 == rhits[2].t
    @test 10 == rhits[3].t
    @test 11 == rhits[4].t
    @test 14 == rhits[5].t
    @test 14 == rhits[6].t
    @test 12 == rhits[7].t
end

# multiplicities
sorted_hits = sort(hits, by=h->h.t)
sort!(sorted_hits, by=h->h.dom_id)
mtps, mtp_ids = NeRCA.count_multiplicities(sorted_hits; tmax=10)
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
f = ROOTFile(ONLINEFILE)
for event in f.online.events
    local hits = NeRCA.combine(event.snapshot_hits, event.triggered_hits)
    @test length(hits) == length(event.snapshot_hits)
    @test length(filter(triggered, hits)) == length(event.triggered_hits)
end

@testset "L1 builder" begin
    f = ROOTFile(datapath("online", "KM3NeT_00000133_JDAQEvent.root"))
    det = Detector(datapath("detx", "KM3NeT_00000133_20221025.detx"))
    l1builder = L1Builder(L1BuilderParameters(25, false))
    @test 33 == length(l1builder(HitL1, det, f.online.events[1].snapshot_hits; combine=false))
    @test 252 == length(l1builder(HitL1, det, f.online.events[2].snapshot_hits; combine=false))
    @test 65 == length(l1builder(HitL1, det, f.online.events[3].snapshot_hits; combine=false))
    @test 10 == length(l1builder(HitL1, det, f.online.events[1].snapshot_hits))
    @test 68 == length(l1builder(HitL1, det, f.online.events[2].snapshot_hits))
    @test 22 == length(l1builder(HitL1, det, f.online.events[3].snapshot_hits))

# Using
# f = NeRCA.ROOTFile("KM3NeT_00000133_JDAQEvent.root")
# det = KM3io.Detector("KM3NeT_00000133_20221025.detx")
# l1builder = NeRCA.L1Builder(NeRCA.L1BuilderParameters(25, false))
#
# This is the reference from Jpp v18
#
# Event #1
# Number of snapshot hits: 590

# L1 hits: 33
# 808443728,89679793.37,2,2 tot=34.292999997735
# 808952499,89678556.434,2,2 tot=41.4350000023842
# 808981913,89681009.963,2,2 tot=36.8369999974966
# 808982006,89680677.946,11,11 tot=52.609999999404
# 808982006,89680680.966,10,10 tot=49.3100000023842
# 808982006,89680681.547,9,9 tot=48.6289999932051
# 808982006,89680681.897,8,8 tot=48.2789999991655
# 808982006,89680681.959,7,7 tot=48.2169999927282
# 808982006,89680682.924,6,6 tot=47.1519999951124
# 808982006,89680682.924,5,5 tot=47.1519999951124
# 808982006,89680682.924,4,4 tot=47.1519999951124
# 808982006,89680682.924,3,3 tot=47.1519999951124
# 808982006,89680686.842,2,2 tot=42.6939999908209
# 808982006,89680696.276,2,2 tot=49.1350000053644
# 808982006,89680718.361,2,2 tot=40.7829999923706
# 809548797,89680703.012,4,4 tot=42.6280000060797
# 809548797,89680710.031,3,3 tot=34.4989999979734
# 809548797,89680718.046,2,2 tot=25.0040000081062
# 816985715,89680915.453,13,13 tot=58.1519999951124
# 816985715,89680917.997,12,12 tot=55.3680000007153
# 816985715,89680920.185,11,11 tot=53
# 816985715,89680919.624,10,10 tot=39.973999992013
# 816985715,89680920.097,9,9 tot=39.5010000020266
# 816985715,89680920.355,8,8 tot=39.2430000007153
# 816985715,89680920.552,7,7 tot=39.0459999889135
# 816985715,89680920.598,6,6 tot=39
# 816985715,89680921.595,5,5 tot=35.0419999957085
# 816985715,89680921.425,4,4 tot=34.2919999957085
# 816985715,89680921.932,3,3 tot=33.625
# 816985715,89680922.372,2,2 tot=33.184999987483
# 817565173,89684035.482,2,2 tot=31.707000002265
# 817589211,89682527.116,2,2 tot=38.402999997139
# 817608176,89684401.156,2,2 tot=42.3050000071526
# L1 combined hits: 10
# 808443728,89679793.37,2,2 tot=34.292999997735
# 808952499,89678556.434,2,2 tot=41.4350000023842
# 808981913,89681009.963,2,2 tot=36.8369999974966
# 808982006,89680677.946,11,11 tot=52.609999999404
# 808982006,89680718.361,2,2 tot=40.7829999923706
# 809548797,89680703.012,4,4 tot=42.6280000060797
# 816985715,89680915.453,13,13 tot=58.1519999951124
# 817565173,89684035.482,2,2 tot=31.707000002265
# 817589211,89682527.116,2,2 tot=38.402999997139
# 817608176,89684401.156,2,2 tot=42.3050000071526
# L2 hits: 6
# 808952499,89678556.434,2,2 tot=41.4350000023842
# 808981913,89681009.963,2,2 tot=36.8369999974966
# 808982006,89680677.946,11,11 tot=52.609999999404
# 809548797,89680703.012,4,4 tot=42.6280000060797
# 816985715,89680915.453,13,13 tot=58.1519999951124
# 817565173,89684035.482,2,2 tot=31.707000002265
#   Match3B on Event #1
# Match3B hits (from L1 combined hits above): road width 200m, tmaxlocal 18.0
# cleaned up hits from L1 combined hits: 9
# 808443728,89679793.37,2,2 tot=34.292999997735
# 808952499,89678556.434,2,2 tot=41.4350000023842
# 808981913,89681009.963,2,2 tot=36.8369999974966
# 808982006,89680677.946,11,11 tot=52.609999999404
# 809548797,89680703.012,4,4 tot=42.6280000060797
# 816985715,89680915.453,13,13 tot=58.1519999951124
# 817565173,89684035.482,2,2 tot=31.707000002265
# 817589211,89682527.116,2,2 tot=38.402999997139
# 817608176,89684401.156,2,2 tot=42.3050000071526
# Match3B hits from L1 combined hits: 5
# 816985715,89680915.453,13,13 tot=58.1519999951124
# 817589211,89682527.116,2,2 tot=38.402999997139
# 808981913,89681009.963,2,2 tot=36.8369999974966
# 808982006,89680677.946,11,11 tot=52.609999999404
# 809548797,89680703.012,4,4 tot=42.6280000060797
# Event #2
# Number of snapshot hits: 1079
# L1 hits: 252
# 806468834,48812232.824,5,5 tot=39.3350000008941
# 806468834,48812233.11,4,4 tot=39.0490000024438
# 806468834,48812233.939,3,3 tot=38.0899999961257
# 806468834,48812234.625,2,2 tot=37.2639999985695
# 808430449,48813184.172,2,2 tot=26
# 808437870,48814613.625,2,2 tot=31.0430000051856
# 808447066,48812517.566,2,2 tot=35.7519999966025
# 808447094,48816096.56,2,2 tot=29
# 808447104,48811799.213,15,15 tot=80
# 808447104,48811803.706,16,16 tot=68.2820000052452
# 808447104,48811805.312,15,15 tot=66.5360000059009
# 808447104,48811805.82,14,14 tot=65.957999996841
# 808447104,48811806.809,13,13 tot=64.8990000039339
# 808447104,48811807.104,12,12 tot=64.6040000021458
# 808447104,48811807.338,11,11 tot=64.3700000047684
# 808447104,48811808.468,10,10 tot=63.1600000038743
# 808447104,48811810.026,9,9 tot=61.4519999995828
# 808447104,48811810.521,8,8 tot=60.8769999966025
# 808447104,48811810.975,7,7 tot=60.4230000078678
# 808447104,48811813.238,6,6 tot=58
# 808447104,48811815.073,5,5 tot=30.9380000010133
# 808447104,48811815.444,4,4 tot=30.5669999942183
# 808447104,48811815.444,3,3 tot=30.5669999942183
# 808447104,48811819.318,2,2 tot=26.8129999935627
# 808447104,48811824.501,2,2 tot=33.9310000017285
# 808447104,48811908.078,2,2 tot=62.3419999927282
# 808447104,48812015.936,2,2 tot=32.9600000008941
# 808451689,48811984.382,14,14 tot=66.5389999970794
# 808451689,48811986.069,13,13 tot=64.7120000049472
# 808451689,48811987.701,12,12 tot=63
# 808451689,48811985.821,11,11 tot=38.4520000070333
# 808451689,48811985.875,11,11 tot=38.3980000019073
# 808451689,48811985.978,10,10 tot=38.2950000017881
# 808451689,48811986.578,9,9 tot=37.554999999702
# 808451689,48811987.06,8,8 tot=37.0729999989271
# 808451689,48811987.865,7,7 tot=36.1180000007153
# 808451689,48811990.114,6,6 tot=33.3890000060201
# 808451689,48811990.486,5,5 tot=33.0170000046492
# 808451689,48811990.316,4,4 tot=32.7470000013709
# 808451689,48811993.036,3,3 tot=29.5469999983907
# 808451689,48811996.943,3,3 tot=46.8559999987483
# 808451689,48812002.601,2,2 tot=40.5179999992251
# 808451689,48812018.439,3,3 tot=29.0890000015497
# 808451689,48812026.208,3,3 tot=49.9499999955297
# 808451689,48812027.088,2,2 tot=48.9699999988079
# 808469556,48811647.105,2,2 tot=43.8259999975562
# 808474243,48812371.267,11,11 tot=56.3379999995232
# 808474243,48812374.306,10,10 tot=52.9490000009537
# 808474243,48812375.218,10,10 tot=51.9469999969006
# 808474243,48812375.617,9,9 tot=51.5480000004172
# 808474243,48812375.649,8,8 tot=51.515999995172
# 808474243,48812377.094,8,8 tot=53.624000005424
# 808474243,48812378.245,7,7 tot=52.3830000013113
# 808474243,48812378.584,6,6 tot=52.0439999997616
# 808474243,48812382.063,5,5 tot=48.1849999949336
# 808474243,48812384.575,4,4 tot=45.3730000033975
# 808474243,48812390.535,4,4 tot=54.8919999971986
# 808474243,48812398.966,4,4 tot=55.1140000000596
# 808474243,48812401.628,4,4 tot=52.1819999963045
# 808474243,48812411.107,3,3 tot=41.5630000010133
# 808474243,48812418.52,2,2 tot=33
# 808474244,48812155.52,2,2 tot=31.8830000013113
# 808493106,48817010.802,2,2 tot=26
# 808947003,48812109.656,2,2 tot=59.1440000012517
# 808947003,48812115,3,3 tot=94.6649999991059
# 808947003,48812141.262,3,3 tot=66.7330000028014
# 808947003,48812141.995,2,2 tot=66
# 808959403,48812114.954,3,3 tot=23.4450000077486
# 808959403,48812122.771,2,2 tot=17.0180000066757
# 808959403,48812191.344,2,2 tot=46.2410000041127
# 808964322,48812146.368,7,7 tot=49.4299999922514
# 808964322,48812147.616,6,6 tot=47.9819999933243
# 808964322,48812147.912,5,5 tot=47.6859999969602
# 808964322,48812148.469,4,4 tot=47.1290000006557
# 808964322,48812148.369,3,3 tot=46.6990000009537
# 808964322,48812162.258,3,3 tot=42.4460000023246
# 808964322,48812165.258,2,2 tot=39.0760000050068
# 808964322,48812184.724,2,2 tot=52.917999997735
# 808964322,48812206.452,2,2 tot=31.7400000020862
# 808964808,48811842.725,6,6 tot=46.9159999936819
# 808964808,48811849.831,5,5 tot=39
# 808964808,48811850.505,4,4 tot=38.1440000012517
# 808964808,48811850.505,3,3 tot=38.1440000012517
# 808964808,48811856.874,2,2 tot=30.5850000008941
# 808964808,48811998.389,2,2 tot=36.8519999980927
# 808964913,48812099.457,3,3 tot=46.1459999978542
# 808964913,48812104.997,2,2 tot=40.33599999547
# 808964925,48812289.832,14,14 tot=55.9689999967813
# 808964925,48812290.063,13,13 tot=55.737999998033
# 808964925,48812290.199,12,12 tot=55.6019999980927
# 808964925,48812290.441,11,11 tot=55.3599999919534
# 808964925,48812292.852,10,10 tot=52.6789999976754
# 808964925,48812293.37,9,9 tot=52.1609999984503
# 808964925,48812293.423,8,8 tot=52.1079999953508
# 808964925,48812293.531,7,7 tot=52
# 808964925,48812303.951,6,6 tot=35.1089999973774
# 808964925,48812304.06,5,5 tot=35
# 808964925,48812304.286,4,4 tot=25.4640000015497
# 808964925,48812304.286,3,3 tot=25.4640000015497
# 808964925,48812310.753,3,3 tot=43.7700000032783
# 808964925,48812311.46,2,2 tot=43.0630000010133
# 808965918,48807378.363,2,2 tot=32.7019999995828
# 808971766,48811867.471,23,23 tot=138.851999990642
# 808971766,48811869.106,22,22 tot=137.186999991536
# 808971766,48811869.552,21,21 tot=136.71099999547
# 808971766,48811869.684,20,20 tot=136.578999996185
# 808971766,48811870.383,19,19 tot=135.84999999404
# 808971766,48811870.548,18,18 tot=135.684999994934
# 808971766,48811870.577,17,17 tot=135.655999995768
# 808971766,48811870.7,16,16 tot=135.532999999821
# 808971766,48811871.145,15,15 tot=135.087999992073
# 808971766,48811871.481,14,14 tot=134.72199999541
# 808971766,48811871.868,13,13 tot=134.334999993443
# 808971766,48811872.175,12,12 tot=133.997999995947
# 808971766,48811872.444,11,11 tot=133.729000002146
# 808971766,48811872.941,10,10 tot=133.231999993324
# 808971766,48811873.173,9,9 tot=133
# 808971766,48811873.14,8,8 tot=128.497000001371
# 808971766,48811873.231,7,7 tot=128.406000003219
# 808971766,48811873.296,6,6 tot=128.341000005603
# 808971766,48811873.637,5,5 tot=128
# 808971766,48811873.167,4,4 tot=100
# 808971766,48811871.472,3,3 tot=27.9159999936819
# 808971766,48811877.278,4,4 tot=48.2820000052452
# 808971766,48811884.241,4,4 tot=40.4390000030398
# 808971766,48811895.162,5,5 tot=46.8229999989271
# 808971766,48811900.68,4,4 tot=40.625
# 808971766,48811906.322,3,3 tot=33.9530000016093
# 808971766,48811910.625,2,2 tot=29
# 808971766,48811940.71,3,3 tot=44.1160000041127
# 808971766,48811952.353,4,4 tot=40.2480000033975
# 808971766,48811955.806,3,3 tot=36.2449999973178
# 808971766,48811963.472,2,2 tot=27.1190000027418
# 808971766,48812013.108,2,2 tot=44.679000005126
# 808972687,48811961.906,2,2 tot=44.9430000036955
# 808973936,48812203.859,3,3 tot=47.2889999970794
# 808973936,48812208.448,2,2 tot=42.1599999964237
# 808976351,48812162.423,2,2 tot=33.3509999960661
# 808978557,48812579.37,2,2 tot=23.4469999969006
# 808979913,48812322.504,4,4 tot=34.4060000032187
# 808979913,48812326.26,3,3 tot=30
# 808979913,48812326.938,2,2 tot=27
# 808984591,48812818.771,2,2 tot=38.0969999954104
# 808984717,48812856.932,2,2 tot=37.4969999939203
# 808984717,48812881.149,2,2 tot=29.4270000010729
# 808987789,48812679.163,2,2 tot=45.9689999967813
# 808998824,48812300.591,2,2 tot=38.3390000015497
# 809538805,48812196.044,2,2 tot=23
# 816964122,48812014.495,4,4 tot=39.0780000016093
# 816964122,48812014.973,3,3 tot=38.4700000062585
# 816964122,48812018.806,2,2 tot=33.8669999986887
# 816964122,48812048.351,2,2 tot=37.4980000033975
# 816964122,48812061.369,2,2 tot=43.6160000041127
# 817318104,48811667.867,4,4 tot=36.7159999981523
# 817318104,48811671.614,3,3 tot=32.3189999982715
# 817318104,48811680.214,2,2 tot=21.8290000036359
# 817333853,48811617.993,2,2 tot=38.7370000034571
# 817333853,48811638.57,3,3 tot=45.4739999994636
# 817333853,48811649.98,2,2 tot=32.2840000092983
# 817333853,48811845.686,3,3 tot=42.5839999914169
# 817333853,48811851.462,2,2 tot=36.0179999917746
# 817565173,48812013.251,19,19 tot=92.3809999972582
# 817565173,48812013.764,18,18 tot=91.8179999962449
# 817565173,48812013.879,17,17 tot=91.7029999941587
# 817565173,48812016.271,17,17 tot=89.2110000029206
# 817565173,48812016.591,16,16 tot=88.8309999927878
# 817565173,48812016.896,15,15 tot=88.5260000005364
# 817565173,48812017.365,14,14 tot=88.0569999963045
# 817565173,48812017.405,13,13 tot=87.9670000001788
# 817565173,48812017.405,12,12 tot=87.9670000001788
# 817565173,48812017.405,11,11 tot=87.9670000001788
# 817565173,48812019.069,10,10 tot=86.2430000007153
# 817565173,48812019.488,9,9 tot=85.7739999964833
# 817565173,48812019.918,8,8 tot=85.3439999967813
# 817565173,48812020.262,7,7 tot=85
# 817565173,48812017.345,6,6 tot=43.554999999702
# 817565173,48812017.9,5,5 tot=43
# 817565173,48812018.588,4,4 tot=38.3150000050664
# 817565173,48812019.727,3,3 tot=37.0360000059009
# 817565173,48812024.803,2,2 tot=31
# 817565173,48812041.129,2,2 tot=45.7279999926686
# 817565173,48812054.637,2,2 tot=43.2520000040531
# 817565173,48812071.509,2,2 tot=48.5459999963641
# 817565173,48812094.775,2,2 tot=37.6000000014901
# 817595425,48811880.557,2,2 tot=37
# 817597172,48809416.939,2,2 tot=34.2390000000596
# 817606442,48811784.386,3,3 tot=39.5059999898076
# 817606442,48811788.322,2,2 tot=35
# 817606442,48811837.352,9,9 tot=44.0509999990463
# 817606442,48811840.245,8,8 tot=40.6880000010133
# 817606442,48811842.094,7,7 tot=38.5789999961853
# 817606442,48811842.935,6,6 tot=37.597999997437
# 817606442,48811843.786,5,5 tot=36.5969999954104
# 817606442,48811844.12,4,4 tot=36.2629999965429
# 817606442,48811844.281,3,3 tot=35.9519999995828
# 817606442,48811848.498,2,2 tot=31.0749999955297
# 817606442,48812035.731,2,2 tot=26
# 817612559,48811899.187,12,12 tot=74.1710000038147
# 817612559,48811900.562,11,11 tot=72.6660000011325
# 817612559,48811901.163,10,10 tot=72.0650000050664
# 817612559,48811902.019,9,9 tot=71.1490000039339
# 817612559,48811902.357,8,8 tot=70.7410000041127
# 817612559,48811903.547,7,7 tot=69.4809999987483
# 817612559,48811903.978,6,6 tot=69.0499999970198
# 817612559,48811904.958,5,5 tot=68
# 817612559,48811902.302,4,4 tot=34.417999997735
# 817612559,48811903.569,3,3 tot=32.820999994874
# 817612559,48811904.981,2,2 tot=31.2489999979734
# 817612559,48811912.76,3,3 tot=51.6380000039935
# 817612559,48811932.822,5,5 tot=37.7669999971986
# 817612559,48811932.822,4,4 tot=37.7669999971986
# 817612559,48811939.438,5,5 tot=48.9490000009537
# 817612559,48811940.167,4,4 tot=48.2199999988079
# 817612559,48811955.589,3,3 tot=30.3980000019073
# 817612559,48811957.397,2,2 tot=28
# 817612559,48811994.511,2,2 tot=38.0370000004768
# 817612559,48812011.458,2,2 tot=25.4839999973774
# 817620192,48811634.474,2,2 tot=29.6620000004768
# 817620193,48812154.601,7,7 tot=49.0200000032783
# 817620193,48812157.01,6,6 tot=46.3109999969602
# 817620193,48812158.894,5,5 tot=44.2170000001788
# 817620193,48812161.426,4,4 tot=41.3350000008941
# 817620193,48812164.381,3,3 tot=38
# 817620193,48812168.276,3,3 tot=52.8189999982715
# 817620193,48812171.517,2,2 tot=49.2980000004172
# 817620193,48812190.435,2,2 tot=37.5949999988079
# 817801282,48811786.386,2,2 tot=26
# 817801282,48811906,4,4 tot=33.0939999967813
# 817801282,48811907.633,3,3 tot=30.9709999933839
# 817801282,48811908.184,2,2 tot=30.4200000017881
# 819006535,48811284.581,17,17 tot=76.1169999986887
# 819006535,48811286.316,16,16 tot=74.2520000040531
# 819006535,48811288.509,15,15 tot=71.8689999952912
# 819006535,48811289.048,14,14 tot=71.3300000056624
# 819006535,48811289.734,14,14 tot=70.5740000009537
# 819006535,48811291.072,13,13 tot=69.1660000085831
# 819006535,48811291.238,12,12 tot=69
# 819006535,48811290.252,11,11 tot=53.4759999960661
# 819006535,48811291.919,10,10 tot=51.6290000006557
# 819006535,48811293.45,9,9 tot=50.0080000013113
# 819006535,48811297.469,9,9 tot=46.320000000298
# 819006535,48811297.724,8,8 tot=46.0649999976158
# 819006535,48811298.453,7,7 tot=45.2360000014305
# 819006535,48811298.826,6,6 tot=44.7529999986291
# 819006535,48811299.503,5,5 tot=44.0759999975562
# 819006535,48811299.903,4,4 tot=43.5659999996424
# 819006535,48811302.378,3,3 tot=40.7309999987483
# 819006535,48811302.699,2,2 tot=40.4099999964237
# 819006535,48811318.109,2,2 tot=40.8339999988675
# 819006535,48811340.353,2,2 tot=32.0720000043511
# 819006535,48811348.405,3,3 tot=47.7170000001788
# 819006535,48811366.353,2,2 tot=26.8790000006557
# L1 combined hits: 68
# 806468834,48812232.824,5,5 tot=39.3350000008941
# 808430449,48813184.172,2,2 tot=26
# 808437870,48814613.625,2,2 tot=31.0430000051856
# 808447066,48812517.566,2,2 tot=35.7519999966025
# 808447094,48816096.56,2,2 tot=29
# 808447104,48811799.213,15,15 tot=80
# 808447104,48811819.318,2,2 tot=26.8129999935627
# 808447104,48811908.078,2,2 tot=62.3419999927282
# 808447104,48812015.936,2,2 tot=32.9600000008941
# 808451689,48811984.382,14,14 tot=66.5389999970794
# 808451689,48812018.439,3,3 tot=29.0890000015497
# 808469556,48811647.105,2,2 tot=43.8259999975562
# 808474243,48812371.267,11,11 tot=56.3379999995232
# 808474243,48812398.966,4,4 tot=55.1140000000596
# 808474244,48812155.52,2,2 tot=31.8830000013113
# 808493106,48817010.802,2,2 tot=26
# 808947003,48812109.656,2,2 tot=59.1440000012517
# 808947003,48812141.262,3,3 tot=66.7330000028014
# 808959403,48812114.954,3,3 tot=23.4450000077486
# 808959403,48812191.344,2,2 tot=46.2410000041127
# 808964322,48812146.368,7,7 tot=49.4299999922514
# 808964322,48812184.724,2,2 tot=52.917999997735
# 808964808,48811842.725,6,6 tot=46.9159999936819
# 808964808,48811998.389,2,2 tot=36.8519999980927
# 808964913,48812099.457,3,3 tot=46.1459999978542
# 808964925,48812289.832,14,14 tot=55.9689999967813
# 808965918,48807378.363,2,2 tot=32.7019999995828
# 808971766,48811867.471,23,23 tot=138.851999990642
# 808971766,48811895.162,5,5 tot=46.8229999989271
# 808971766,48811940.71,3,3 tot=44.1160000041127
# 808971766,48811963.472,2,2 tot=27.1190000027418
# 808971766,48812013.108,2,2 tot=44.679000005126
# 808972687,48811961.906,2,2 tot=44.9430000036955
# 808973936,48812203.859,3,3 tot=47.2889999970794
# 808976351,48812162.423,2,2 tot=33.3509999960661
# 808978557,48812579.37,2,2 tot=23.4469999969006
# 808979913,48812322.504,4,4 tot=34.4060000032187
# 808984591,48812818.771,2,2 tot=38.0969999954104
# 808984717,48812856.932,2,2 tot=37.4969999939203
# 808987789,48812679.163,2,2 tot=45.9689999967813
# 808998824,48812300.591,2,2 tot=38.3390000015497
# 809538805,48812196.044,2,2 tot=23
# 816964122,48812014.495,4,4 tot=39.0780000016093
# 816964122,48812048.351,2,2 tot=37.4980000033975
# 817318104,48811667.867,4,4 tot=36.7159999981523
# 817333853,48811617.993,2,2 tot=38.7370000034571
# 817333853,48811649.98,2,2 tot=32.2840000092983
# 817333853,48811845.686,3,3 tot=42.5839999914169
# 817565173,48812013.251,19,19 tot=92.3809999972582
# 817565173,48812041.129,2,2 tot=45.7279999926686
# 817565173,48812071.509,2,2 tot=48.5459999963641
# 817595425,48811880.557,2,2 tot=37
# 817597172,48809416.939,2,2 tot=34.2390000000596
# 817606442,48811784.386,3,3 tot=39.5059999898076
# 817606442,48811837.352,9,9 tot=44.0509999990463
# 817606442,48812035.731,2,2 tot=26
# 817612559,48811899.187,12,12 tot=74.1710000038147
# 817612559,48811932.822,5,5 tot=37.7669999971986
# 817612559,48811957.397,2,2 tot=28
# 817612559,48811994.511,2,2 tot=38.0370000004768
# 817620192,48811634.474,2,2 tot=29.6620000004768
# 817620193,48812154.601,7,7 tot=49.0200000032783
# 817620193,48812190.435,2,2 tot=37.5949999988079
# 817801282,48811786.386,2,2 tot=26
# 817801282,48811906,4,4 tot=33.0939999967813
# 819006535,48811284.581,17,17 tot=76.1169999986887
# 819006535,48811318.109,2,2 tot=40.8339999988675
# 819006535,48811348.405,3,3 tot=47.7170000001788
# L2 hits: 55
# 806468834,48812232.824,5,5 tot=39.3350000008941
# 808437870,48814613.625,2,2 tot=31.0430000051856
# 808447066,48812517.566,2,2 tot=35.7519999966025
# 808447094,48816096.56,2,2 tot=29
# 808447104,48811799.213,15,15 tot=80
# 808447104,48811819.318,2,2 tot=26.8129999935627
# 808447104,48811908.078,2,2 tot=62.3419999927282
# 808451689,48811984.382,14,14 tot=66.5389999970794
# 808451689,48812018.439,3,3 tot=29.0890000015497
# 808474243,48812371.267,11,11 tot=56.3379999995232
# 808474243,48812398.966,4,4 tot=55.1140000000596
# 808493106,48817010.802,2,2 tot=26
# 808947003,48812109.656,2,2 tot=59.1440000012517
# 808947003,48812141.262,3,3 tot=66.7330000028014
# 808959403,48812114.954,3,3 tot=23.4450000077486
# 808964322,48812146.368,7,7 tot=49.4299999922514
# 808964322,48812184.724,2,2 tot=52.917999997735
# 808964808,48811842.725,6,6 tot=46.9159999936819
# 808964808,48811998.389,2,2 tot=36.8519999980927
# 808964913,48812099.457,3,3 tot=46.1459999978542
# 808964925,48812289.832,14,14 tot=55.9689999967813
# 808971766,48811867.471,23,23 tot=138.851999990642
# 808971766,48811895.162,5,5 tot=46.8229999989271
# 808971766,48811940.71,3,3 tot=44.1160000041127
# 808973936,48812203.859,3,3 tot=47.2889999970794
# 808976351,48812162.423,2,2 tot=33.3509999960661
# 808978557,48812579.37,2,2 tot=23.4469999969006
# 808979913,48812322.504,4,4 tot=34.4060000032187
# 808984591,48812818.771,2,2 tot=38.0969999954104
# 808987789,48812679.163,2,2 tot=45.9689999967813
# 808998824,48812300.591,2,2 tot=38.3390000015497
# 809538805,48812196.044,2,2 tot=23
# 816964122,48812014.495,4,4 tot=39.0780000016093
# 816964122,48812048.351,2,2 tot=37.4980000033975
# 817318104,48811667.867,4,4 tot=36.7159999981523
# 817333853,48811617.993,2,2 tot=38.7370000034571
# 817333853,48811649.98,2,2 tot=32.2840000092983
# 817333853,48811845.686,3,3 tot=42.5839999914169
# 817565173,48812013.251,19,19 tot=92.3809999972582
# 817565173,48812041.129,2,2 tot=45.7279999926686
# 817565173,48812071.509,2,2 tot=48.5459999963641
# 817595425,48811880.557,2,2 tot=37
# 817597172,48809416.939,2,2 tot=34.2390000000596
# 817606442,48811784.386,3,3 tot=39.5059999898076
# 817606442,48811837.352,9,9 tot=44.0509999990463
# 817606442,48812035.731,2,2 tot=26
# 817612559,48811899.187,12,12 tot=74.1710000038147
# 817612559,48811932.822,5,5 tot=37.7669999971986
# 817612559,48811957.397,2,2 tot=28
# 817620192,48811634.474,2,2 tot=29.6620000004768
# 817620193,48812154.601,7,7 tot=49.0200000032783
# 817801282,48811786.386,2,2 tot=26
# 817801282,48811906,4,4 tot=33.0939999967813
# 819006535,48811284.581,17,17 tot=76.1169999986887
# 819006535,48811318.109,2,2 tot=40.8339999988675
# Event #3
# Number of snapshot hits: 697
# L1 hits: 65
# 806459597,81498216.734,2,2 tot=37.5489999949932
# 808443728,81497475.73,2,2 tot=44.777999997139
# 808474578,81496601.079,2,2 tot=26.4800000041723
# 808474578,81496702.355,2,2 tot=47.9050000011921
# 808961256,81497949.399,2,2 tot=26
# 808977112,81497699.319,2,2 tot=31.1409999877214
# 809003775,81495714.643,2,2 tot=34.7189999967813
# 816919154,81497186.185,5,5 tot=52.6070000082254
# 816919154,81497186.185,4,4 tot=52.6070000082254
# 816919154,81497188.319,3,3 tot=50.2930000126362
# 816919154,81497188.319,2,2 tot=50.2930000126362
# 816919154,81497229.819,2,2 tot=26.4230000078678
# 816930809,81497207.622,6,6 tot=59.5259999930859
# 816930809,81497207.63,5,5 tot=59.5179999917746
# 816930809,81497207.63,4,4 tot=59.5179999917746
# 816930809,81497208.976,3,3 tot=58.0919999927282
# 816930809,81497209.068,2,2 tot=58
# 816930809,81497214.302,2,2 tot=46.2999999970198
# 816930809,81497238.042,2,2 tot=41.9320000112057
# 816985715,81496569.392,12,12 tot=63.3199999928474
# 816985715,81496569.68,11,11 tot=62.9619999974966
# 816985715,81496569.68,10,10 tot=62.9619999974966
# 816985715,81496572.529,9,9 tot=59.8729999959469
# 816985715,81496574.408,9,9 tot=57.8340000063181
# 816985715,81496574.712,8,8 tot=57.5300000011921
# 816985715,81496574.76,7,7 tot=57.4819999933243
# 816985715,81496575.159,6,6 tot=57.0829999893904
# 816985715,81496575.242,5,5 tot=57
# 816985715,81496574.521,4,4 tot=42.1120000034571
# 816985715,81496577.263,3,3 tot=39
# 816985715,81496580.609,2,2 tot=35
# 816985715,81496600.365,2,2 tot=33.397000014782
# 817287557,81496364.93,14,14 tot=95.4579999893904
# 817287557,81496365.267,13,13 tot=95.1210000067949
# 817287557,81496366.874,12,12 tot=93.41400000453
# 817287557,81496367.062,11,11 tot=93.2260000109673
# 817287557,81496367.109,10,10 tot=93.1789999902248
# 817287557,81496367.527,9,9 tot=92.7110000103712
# 817287557,81496368.03,8,8 tot=92.2080000042915
# 817287557,81496368.03,7,7 tot=92.2080000042915
# 817287557,81496365.25,6,6 tot=51
# 817287557,81496365.25,5,5 tot=51
# 817287557,81496365.25,4,4 tot=51
# 817287557,81496365.25,3,3 tot=51
# 817287557,81496366.49,2,2 tot=37
# 817287557,81496397.779,2,2 tot=41.7910000085831
# 817318910,81496424.004,6,6 tot=45
# 817318910,81496424.657,5,5 tot=43.5509999990463
# 817318910,81496427.94,5,5 tot=57.1229999959469
# 817318910,81496428.593,4,4 tot=56.390000000596
# 817318910,81496432.608,3,3 tot=51.9350000023842
# 817318910,81496439.321,2,2 tot=44.5219999998808
# 817318910,81496487.723,1,1 tot=4
# 817344240,81496633.579,4,4 tot=52
# 817344240,81496633.579,3,3 tot=52
# 817344240,81496635.09,2,2 tot=44
# 817802243,81496794.389,2,2 tot=44.847000002861
# 819005766,81500561.171,2,2 tot=44.3879999965429
# 819006150,81496715.615,2,2 tot=41.5539999902248
# 819047388,81497219.983,6,6 tot=51.0130000114441
# 819047388,81497220.643,5,5 tot=50.2630000114441
# 819047388,81497223.785,4,4 tot=46.7210000008345
# 819047388,81497232.749,4,4 tot=40.7680000066757
# 819047388,81497235.519,3,3 tot=37.5980000048876
# 819047388,81497239.626,2,2 tot=32.6909999996424
# L1 combined hits: 22
# 806459597,81498216.734,2,2 tot=37.5489999949932
# 808443728,81497475.73,2,2 tot=44.777999997139
# 808474578,81496601.079,2,2 tot=26.4800000041723
# 808474578,81496702.355,2,2 tot=47.9050000011921
# 808961256,81497949.399,2,2 tot=26
# 808977112,81497699.319,2,2 tot=31.1409999877214
# 809003775,81495714.643,2,2 tot=34.7189999967813
# 816919154,81497186.185,5,5 tot=52.6070000082254
# 816919154,81497229.819,2,2 tot=26.4230000078678
# 816930809,81497207.622,6,6 tot=59.5259999930859
# 816930809,81497238.042,2,2 tot=41.9320000112057
# 816985715,81496569.392,12,12 tot=63.3199999928474
# 816985715,81496600.365,2,2 tot=33.397000014782
# 817287557,81496364.93,14,14 tot=95.4579999893904
# 817287557,81496397.779,2,2 tot=41.7910000085831
# 817318910,81496424.004,6,6 tot=45
# 817318910,81496487.723,1,1 tot=4
# 817344240,81496633.579,4,4 tot=52
# 817802243,81496794.389,2,2 tot=44.847000002861
# 819005766,81500561.171,2,2 tot=44.3879999965429
# 819006150,81496715.615,2,2 tot=41.5539999902248
# 819047388,81497219.983,6,6 tot=51.0130000114441
# L2 hits: 17
# 808443728,81497475.73,2,2 tot=44.777999997139
# 808474578,81496601.079,2,2 tot=26.4800000041723
# 808474578,81496702.355,2,2 tot=47.9050000011921
# 808977112,81497699.319,2,2 tot=31.1409999877214
# 809003775,81495714.643,2,2 tot=34.7189999967813
# 816919154,81497186.185,5,5 tot=52.6070000082254
# 816919154,81497229.819,2,2 tot=26.4230000078678
# 816930809,81497207.622,6,6 tot=59.5259999930859
# 816930809,81497238.042,2,2 tot=41.9320000112057
# 816985715,81496569.392,12,12 tot=63.3199999928474
# 817287557,81496364.93,14,14 tot=95.4579999893904
# 817287557,81496397.779,2,2 tot=41.7910000085831
# 817318910,81496424.004,6,6 tot=45
# 817344240,81496633.579,4,4 tot=52
# 819005766,81500561.171,2,2 tot=44.3879999965429
# 819006150,81496715.615,2,2 tot=41.5539999902248
# 819047388,81497219.983,6,6 tot=51.0130000114441
end
