using KM3NeT
using Test


# angle_between()
@test 0 == KM3NeT.angle_between([1,0,0], [1,0,0])
@test π/2 ≈ KM3NeT.angle_between([1,0,0], [0,1,0])
@test π/2 ≈ KM3NeT.angle_between([1,0,0], [0,0,1])
@test π ≈ KM3NeT.angle_between([1,0,0], [-1,0,0])

# make_cherenkov_calculator()
# ccalc = KM3NeT.make_cherenkov_calculator(KM3NeT.Track([0,0,-1], [0,0,1], 0))
# @test 3.33564095 =≈ ccalc([0, 0, 0])
