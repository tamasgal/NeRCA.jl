using KM3NeT
using Test


# make_cherenkov_calculator()
ccalc = KM3NeT.make_cherenkov_calculator([0,0,1], [0,0,-1])
@test 3.33564095 ≈ ccalc([0, 0, 0])
