using NeRCA
using KM3io
using Test


@testset "Single DU" begin
    f = ROOTFile(joinpath(@__DIR__, "data", "mcv6.1.mupage_10G.sirene.jterbr00007209.2548_0-499_JDAQEvent.root"))
    det = Detector(joinpath(@__DIR__, "data", "KM3NeT_00000044_00007209.v6.0_PMTeff_merge8.K40.detx"))
    sparams = NeRCA.SingleDURecoParams(floor_distance=floordist(det))

    # event_id => (brightest, reco_dz)
    results = Dict(
        1 => (3, -0.736740691005346),
        2 => (3, -0.9500944406773869),
        3 => (5, -0.9290231912579134),
        23 => (4, -0.6062414827877373),
        500 => (5, -0.99119883292077),
    )

    for event_id in keys(results)
        event = f.online.events[event_id]

        hits = event.snapshot_hits
        thits = event.triggered_hits

        chits = calibrate(det, NeRCA.combine(hits, thits))

        brightest_du = most_frequent(h->h.du, triggered(chits))
        @test results[event_id][1] == brightest_du

        du_hits = filter(h->h.du == brightest_du, chits)

        fit = NeRCA.single_du_fit(du_hits, sparams)
        # @test results[event_id][2] ≈ fit.sdp.dz
    end

end
