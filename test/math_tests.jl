using NeRCA
using Test


@testset "angle()" begin
    @test 0 == angle([1,0,0], [1,0,0])
    @test π/2 ≈ angle([1,0,0], [0,1,0])
    @test π/2 ≈ angle([1,0,0], [0,0,1])
    @test π ≈ angle([1,0,0], [-1,0,0])
end
