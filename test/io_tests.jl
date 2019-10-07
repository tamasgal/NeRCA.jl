using NeRCA
using Test


const DETX = joinpath(@__DIR__, "data", "detx_v3.detx")


@testset "calibration" begin
    calib = NeRCA.read_calibration(DETX)

    @test 23 == calib.det_id
    @test 6 == length(values(calib.pos))
    @test 6 == length(values(calib.dir))
    @test 6 == length(values(calib.t0))
    @test 6 == length(values(calib.du))
    @test 6 == length(values(calib.floor))
    @test 18 == length(values(calib.omkeys))
    @test 23.9 ≈ calib.max_z
    @test 2 == calib.n_dus
    @test 3.6 ≈ calib.pos[3][2].z
end
