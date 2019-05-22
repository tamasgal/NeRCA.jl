using KM3NeT
using Test


# angle_between()
@test 0 == KM3NeT.angle_between([1,0,0], [1,0,0])
@test π/2 ≈ KM3NeT.angle_between([1,0,0], [0,1,0])
@test π/2 ≈ KM3NeT.angle_between([1,0,0], [0,0,1])
@test π ≈ KM3NeT.angle_between([1,0,0], [-1,0,0])

# azimuth()
@test 0 == KM3NeT.azimuth(Direction(0,1,0))
@test π/2 ≈ KM3NeT.azimuth(Direction(1,0,0))
@test π ≈ KM3NeT.azimuth(Direction(0,-1,0))
@test (3/2)π ≈ KM3NeT.azimuth(Direction(-1,0,0))
