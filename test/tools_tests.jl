using NeRCA
using Test


@testset "most_frequent()" begin
    a = [1, 1, 2, 3, 1, 5]
    @test 1 == most_frequent(x->x, a)

    a = [[1, 2], [1, 2, 3], [1, 2], [1, 2], [1], [1]]
    @test 3 == most_frequent(sum, a)
end
