using KM3NeT
using Test


# make_cherenkov_calculator()
ccalc = KM3NeT.make_cherenkov_calculator(KM3NeT.Track([0,0,-1], [0,0,1], 0))
@test 3.33564095 ≈ ccalc([0, 0, 0])
