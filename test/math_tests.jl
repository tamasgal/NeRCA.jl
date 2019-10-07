using NeRCA
using Test


@testset "angle_between()" begin
    @test 0 == NeRCA.angle_between([1,0,0], [1,0,0])
    @test π/2 ≈ NeRCA.angle_between([1,0,0], [0,1,0])
    @test π/2 ≈ NeRCA.angle_between([1,0,0], [0,0,1])
    @test π ≈ NeRCA.angle_between([1,0,0], [-1,0,0])
end

@testset "azimuth()" begin
    @test π/2 == NeRCA.azimuth(Direction(0,1,0))
    @test 0 ≈ NeRCA.azimuth(Direction(1,0,0))
    @test -π/2 ≈ NeRCA.azimuth(Direction(0,-1,0))
    @test π ≈ NeRCA.azimuth(Direction(-1,0,0))
end

@testset "zenith()" begin
    @test 0.0 == NeRCA.zenith(Direction(0,0,-1))
    @test π/2 == NeRCA.zenith(Direction(0,1,0))
    @test π/2 == NeRCA.zenith(Direction(1,1,0))
    @test π ≈ NeRCA.zenith(Direction(0,0,1))
end
