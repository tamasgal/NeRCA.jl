#!/usr/bin/env julia

if length(ARGS) < 2
    println("Usage: ./dumandfit.jl DETECTORFILE ONLINEFILE")
    exit(1)
end

using NeRCA
using KM3io
using ProgressMeter
# using FHist

struct DumandFit
    x::Float32
    y::Float32
    z::Float32
    dx::Float32
    dy::Float32
    dz::Float32
end

struct RecoPerformance
    walltime::Float32
    angular_error::Float32
    mc_zenith::Float32
    reco_zenith::Float32
    mc_nu_energy::Float32
    n_hits_used::Int32
end

function main()
    det = Detector(ARGS[1])
    fname = ARGS[2]
    f = KM3io.ROOTFile(fname)

    outfile = H5File("test.h5", "w")
    reco_dset = create_dataset(outfile, "reco/dumand", DumandFit)
    reco_performance_dset = create_dataset(outfile, "reco_performance/dumand", RecoPerformance)


    @showprogress 1 for event ∈ f.online.events
        thits = calibrate(det, event.triggered_hits)

        walltime = @elapsed track = dumandfit(thits)

        mc_event = f.offline[event.header.trigger_counter + 1]  # one-based indexing in Julia

        ν = first(mc_event.mc_trks)  # primary neutrino by definition the first MC track
        ξ = angle(track.dir, ν.dir)

        rp = RecoPerformance(walltime, ξ, zenith(ν.dir), zenith(track.dir), ν.E, length(thits))

        push!(reco_dset, DumandFit(track.pos..., track.dir...))
        push!(reco_performance_dset, rp)
    end
    close(outfile)
    close(f)
end

main()
