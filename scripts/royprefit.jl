#!/usr/bin/env julia

if length(ARGS) < 2
    println("Usage: ./royprefit.jl DETX HDF5FILE OUTFILE")
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

    outf = open(outfile, "w")
    write(outf, "group_id,dx,dy,dz,x,y,z,t0\n")

    @showprogress 1 for event in KM3NeT.EventReader(filename, detx)
        hits = calibrate(event.hits, event.calib)
        sort!(hits, by=h->h.t)
        sort!(hits, by=h->h.dom_id)
        KM3NeT.count_multiplicities!(hits)

        shits = filter(h->h.multiplicity.count >= 4, hits)
        doms = unique(map(h->h.dom_id, shits))
        track = KM3NeT.prefit(shits)
        write(outf, "$(event.info.group_id),$(track.dir.x),$(track.dir.y),$(track.dir.z),$(track.pos.x),$(track.pos.y),$(track.pos.z),$(track.time)\n")
    end

    close(outf)
end

main()
