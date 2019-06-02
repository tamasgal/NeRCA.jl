#!/usr/bin/env julia

if length(ARGS) < 2
    println("Usage: ./hit_stats.jl DETX HDF5FILE OUTFILE")
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

    events = KM3NeT.EventReader(filename, detx)

    fobj = h5open(filename, "r")

    outf = open(outfile, "w")
    write(outf, "time_residual,n_muons,tot,multiplicity,triggered,distance,bundle_energy,track_pmt_angle\n")

    @showprogress 1 for event in events
        n_muons = length(event.mc_tracks)
        bundle_energy = sum(map(m->m.E, event.mc_tracks))
        hits = calibrate(event.hits, event.calib)
        muon = KM3NeT.Track(event.mc_tracks[1])
        KM3NeT.count_multiplicities!(hits)
        ccalc = KM3NeT.make_cherenkov_calculator(muon)
        mc_time = KM3NeT.make_mc_time_converter(event.info)
        for hit in hits
            Δt = hit.t - mc_time(ccalc(hit.pos))
            d = KM3NeT.pld3(hit.pos, muon.pos, muon.dir)
            η = KM3NeT.angle_between(hit.dir, muon.dir)
            write(outf, "$(Δt),$(n_muons),$(hit.tot),$(hit.multiplicity.count),$(hit.triggered ? 1 : 0),$(d),$(bundle_energy),$(η)\n")
        end
    end

    close(outf)
    close(fobj)
end

main()
