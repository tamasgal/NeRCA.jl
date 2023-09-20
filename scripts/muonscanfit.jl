#!/usr/bin/env julia
using NeRCA
using KM3io
using ProgressMeter

function main()
    fname = ARGS[1]
    f = NeRCA.ROOTFile(fname)
    det = KM3io.Detector("KM3NeT_00000133_20221025.detx")

    msfparams = MuonScanfitParameters(;tmaxlocal=18.0, roadwidth=200.0)
    msfit = MuonScanfit(msfparams, det)

    outfile = open(basename(fname) * ".txt", "w")

    # event 1 (trigger_counter 4)
    #  Direction(0.041586942949292174, -0.029550312319766795, 0.9986978047527373)
    # event 2 (trigger_counter 5)
    #  Direction(0.7284017634046146, -0.28083929562338567, 0.6249481267295778)
    # event 3 (trigger_counter 17)
    # Direction(-0.1129797622905296, -0.484448542891068, 0.8674936210736632)

    # for event_id ∈ [3, 1, 2, 4, 5]  # Jpp readout order, sorted by frame_index
    #     # println("\n\n========================================")
    #     event = f.online.events[event_id]
    #     # @show event_id event.header.frame_index
    #     muons = msfit(event.snapshot_hits)
    #     # @show muons
    # end

    n_events = 0
    n_muons = 0
    n_failed = 0
    @showprogress 1 for (event_id, event) ∈ enumerate(f.online.events)
        muons = msfit(event.snapshot_hits)
        n_events += 1
        n_muons += length(muons)
        if length(muons) == 0
            n_failed += 1
            continue
        end
        best_muon = first(muons)
        neutrino = first(f.offline[event.header.trigger_counter + 1].mc_trks)
        angular_error = angle(best_muon.dir, neutrino.dir)
        write(outfile, "$angular_error\n")
        flush(outfile)
    end
    close(outfile)

    println("Total number of events processed: $n_events")
    println("Number of muon candidates: $n_muons")
    println("Failed reconstructions: $n_failed")
end

main()
