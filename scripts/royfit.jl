#!/usr/bin/env julia

if length(ARGS) < 2
    println("Usage: ./royfit.jl DETX HDF5FILE OUTFILE")
    exit(1)
end


using KM3NeT
using HDF5
using ProgressMeter


function main()
    detx, filename, outfile = ARGS
    if isfile(outfile)
        println("Output file already exists, aborting...")
        exit(1)
    end

    calib = KM3NeT.read_calibration(detx)
    events = KM3NeT.read_compound(filename, "/event_info", KM3NeT.MCEventInfo);
    fobj = h5open(filename, "r")

    outf = open(outfile, "w")
    write(outf, "group_id,d,t,z,dz,phi,t0\n")

    @showprogress 1 for event in events
        hits = calibrate(KM3NeT.read_hits(fobj, event.group_id), calib)
        triggered_hits = filter(h -> h.triggered, hits);
        brightest_du = KM3NeT.most_frequent(h -> h.du, triggered_hits)
        du_hits = filter(h -> h.du == brightest_du, hits)
        fit = KM3NeT.reco(du_hits)
        write(outf, "$(event.group_id),$(fit.sdp.d),$(fit.sdp.t),$(fit.sdp.z),$(fit.sdp.dz),$(fit.sdp.ϕ),$(fit.sdp.t₀)\n")
    end

    close(outf)
    close(fobj)
end

main()
