using NeRCA
using Test


@testset "angle()" begin
    @test 0 == angle([1,0,0], [1,0,0])
    @test π/2 ≈ angle([1,0,0], [0,1,0])
    @test π/2 ≈ angle([1,0,0], [0,0,1])
    @test π ≈ angle([1,0,0], [-1,0,0])

    m1 = NeRCA.MuonScanfitCandidate(Position(1, 2, 3), Direction(1.0, 0.0, 0.0), 0, 0, 0, 0)
    m2 = NeRCA.MuonScanfitCandidate(Position(1, 2, 3), Direction(0.0, 1.0, 0.0), 0, 0, 0, 0)
    @test π/2 ≈ angle(m1, m2)
end


@testset "spread()" begin
    directions = [
        Direction(1.0, 0.0, 0.0),
        Direction(1.0, 0.0, 0.0),
    ]
    @test 0 ≈ spread(directions)
    directions = [
        Direction(1.0, 0.0, 0.0),
        Direction(0.0, 1.0, 0.0),
    ]
    @test π/2 ≈ spread(directions)
    directions = [
        Direction(1.0, 0.0, 0.0),
        Direction(-1.0, 0.0, 0.0),
    ]
    @test π ≈ spread(directions)
end
