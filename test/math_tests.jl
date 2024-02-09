using NeRCA
using Test


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
