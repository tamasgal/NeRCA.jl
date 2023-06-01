#!/usr/bin/env julia

if length(ARGS) < 1
    println("Usage: ./dumandfit.jl OFFLINEFILE")
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
    n_dus::Int32
    angular_error::Float32
    mc_zenith::Float32
    reco_zenith::Float32
    mc_nu_energy::Float32
    n_hits_used::Int32
end

function main()
    fname = ARGS[1]
    f = KM3io.ROOTFile(fname)
    outfile = H5File("test.h5", "w")
    reco_dset = create_dataset(outfile, "reco/dumand", DumandFit)
    reco_performance_dset = create_dataset(outfile, "reco_performance/dumand", RecoPerformance)

    @showprogress 1 for event ∈ f.offline
        thits = filter(triggered, event.hits)

        walltime = @elapsed track = dumandfit(thits)

        ν = first(event.mc_trks)  # primary neutrino by definition the first MC track
        ξ = angle(track.dir, ν.dir)
        n_dus = length(unique([h.dom_id for h ∈ thits]))

        rp = RecoPerformance(walltime, n_dus, ξ, zenith(ν.dir), zenith(track.dir), ν.E, length(thits))

        push!(reco_dset, DumandFit(track.pos..., track.dir...))
        push!(reco_performance_dset, rp)
    end
    close(outfile)
end

main()
