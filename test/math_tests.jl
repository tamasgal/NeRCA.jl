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

@testset "lld3()" begin
    @test 5.0 ≈ NeRCA.lld3(
        Position(0.0, 0.0, 0.0),
        Direction(0.0, 0.0, 1.0),
        Position(5.0, 0.0, 0.0),
        Direction(0.0, 1.0, 0.0),
    )

    @test 4.108795387706403 ≈ NeRCA.lld3(
        Position(1.0, 2.0, 3.0),
        Direction(0.6, 0.3, 1.0),
        Position(-1.0, 5.0, 6.0),
        Direction(0.8, 1.0, 0.3),
    )

    # parallel lines
    @test 4.0 ≈ NeRCA.lld3(
        Position(0.0, 1.0, 0.0),
        Direction(1.0, 0.0, 0.0),
        Position(0.0, 5.0, 0.0),
        Direction(1.0, 0.0, 0.0),
    )
end
