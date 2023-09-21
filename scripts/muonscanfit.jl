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
