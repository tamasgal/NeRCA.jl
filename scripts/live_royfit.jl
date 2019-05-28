#!/usr/bin/env julia
using KM3NeT
using Plots
GR.inline("png")

if length(ARGS) < 1
    println("Usage: ./live_royfit.jl DETX")
    exit(1)
end

const calib = KM3NeT.read_calibration(ARGS[1])


function main()
    println("Starting live ROyFit")
    for message in CHClient(ip"127.0.0.1", 55530, ["IO_EVT"])
        event = KM3NeT.read_io(IOBuffer(message.data), KM3NeT.DAQEvent)

        track = reconstruct(event)
        if track == nothing
            println("Skipping")
            continue
        end
        println(track)
        if rand() > 0.94
            println("Plotting...")
            snapshot_hits = calibrate(event.snapshot_hits, calib)
            triggered_hits = calibrate(event.triggered_hits, calib)
            plot(snapshot_hits)
            plot!(triggered_hits, track)
            savefig("ztplot.png")
        end
    end
end


function reconstruct(event)
    triggered_hits = calibrate(event.triggered_hits, calib)
    n_dus = length(unique(h->h.du, triggered_hits))
    if n_dus < 2
        return nothing
    end
    n_doms = length(unique(h->h.dom_id, triggered_hits))
    if n_doms < 4
        return nothing
    end
    # snapshot_hits = calibrate(event.snapshot_hits, calib)

    return KM3NeT.prefit(triggered_hits)
end


main()
