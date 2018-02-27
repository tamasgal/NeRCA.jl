using KM3NeT
using Base.Test


hits = [Hit(1, 8, 100, 20, false),
        Hit(1, 9, 101, 21, true),
        Hit(2, 8, 112, 22, true),
        Hit(2, 9, 114, 23, false),
        Hit(1, 8, 134, 24, true),
        Hit(2, 9, 156, 25, false),
        Hit(2, 8, 133, 26, true),
        Hit(0, 8, 145, 26, false)]


calib = Calibration(
    # detector ID
    1,
    # PMT positions
    Dict{Int32, Vector{KM3NeT.Position}}(
        8 => [(1,2,3), (4,5,6), (7,8,9)],
        9 => [(10,11,12), (13,14,15),(16,17,18)]
    ),
    # PMT directions
    Dict{Int32, Vector{KM3NeT.Direction}}(
        8 => [(1,2,3), (4,5,6), (7,8,9)],
        9 => [(10,11,12), (13,14,15),(16,17,18)]
    ),
    # t0
    Dict{Int32, Vector{Int32}}(
        8 => [1000,2000,3000],
        9 => [4000,5000,6000]
    ),
    # DU
    Dict{Int32, UInt8}(8 => 42, 9 => 23),
    # floor
    Dict{Int32, UInt8}(8 => 18, 9 => 5),
)


chits = calibrate(hits, calib)
@test 42 == chits[1].du
@test KM3NeT.Position(4,5,6) == chits[1].pos
@test KM3NeT.Direction(16,17,18) == chits[6].dir
@test 5000 == chits[2].t0
@test 18 == chits[3].floor
