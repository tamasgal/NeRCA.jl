using NeRCA
using Test


@testset "most_frequent()" begin
    a = [1, 1, 2, 3, 1, 5]
    @test 1 == most_frequent(x->x, a)

    a = [[1, 2], [1, 2, 3], [1, 2], [1, 2], [1], [1]]
    @test 3 == most_frequent(sum, a)
end


@testset "nthbitset()" begin
    @test NeRCA.nthbitset(2, 12)
    for n ∈ [1, 3, 5]
        @test NeRCA.nthbitset(n, 42)
    end
    for n ∈ [1, 5, 7, 10, 14, 17, 18, 19, 20, 25, 26, 29, 31, 32, 33, 35, 38, 40, 41, 43, 44, 47, 49, 50, 52, 53, 55, 56]
        @test NeRCA.nthbitset(n, 123456789011121314)
    end
end
