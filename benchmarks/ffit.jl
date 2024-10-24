#!/usr/bin/env julia
println("Loading libraries...")
using LinearAlgebra
using NeRCA
using KM3io
using ProgressMeter
using KM3DB

if Threads.nthreads() > 1
  println("Running with $(Threads.nthreads()) Julia threads, reducing BLAS threads to 1 to avoid thread fight.")
  BLAS.set_num_threads(1)
end


struct FFitResult
    det_id::Int
    run::Int
    frame_index::Int
    trigger_counter::Int
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    dir_x::Float64
    dir_y::Float64
    dir_z::Float64
    azimuth::Float64
    zenith::Float64
    t::Float64
    Q::Float64
    S1::Float64
    S2::Float64
    NDF::Int
    walltime::Float64
end


struct MCInfo
    det_id::Int
    run::Int
    frame_index::Int
    trigger_counter::Int
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    dir_x::Float64
    dir_y::Float64
    dir_z::Float64
    azimuth::Float64
    zenith::Float64
    t::Float64
    energy::Float64
    type::Int
    bundle_size::Int
end


function main()
    fname = joinpath(@__DIR__, "data/mcv8.1.mupage_tuned_100G.sirene.jterbr00013288.1.root")
    detfname = joinpath(@__DIR__, "data/KM3NeT_00000133_20221025.detx")
    outfname = "fooooo.h5"

    println("Opening $fname")
    f = ROOTFile(fname)

    println("Loading detector description from $(detfname)")
    det = Detector(detfname)

    msfparams = MuonScanfitParameters(;tmaxlocal=18.0, roadwidth=200.0, nfits=12, α₁=7.0, α₂=0.3, θ=4.0)

    outfile = H5File(outfname, "w")
    reco_dset = create_dataset(outfile, "reco/ffit", FFitResult)
    addmeta(reco_dset, msfparams)

    if f.offline !== nothing
        println("Found MC info => also creating MC histograms")
        mc_dset = create_dataset(outfile, "mc", MCInfo)
    end

    n = Threads.Atomic{Int}(0)
    n_muons = Threads.Atomic{Int}(0)
    n_failed = Threads.Atomic{Int}(0)

    progressbar = Progress(length(f.online.events); desc="Reconstructing: ", dt=0.2, showspeed=true)
    Threads.@threads for event in f.online.events
        next!(progressbar)
        msfit = MuonScanfit(msfparams, det)

        walltime = @elapsed muons = msfit(event.snapshot_hits)

        Threads.atomic_add!(n, 1)
        Threads.atomic_add!(n_muons, muons |> length)
        if length(muons) == 0
          Threads.atomic_add!(n_failed, 1)
          continue
        end

        if f.offline !== nothing
            mc_event = f.offline[event.header.trigger_counter + 1]
            primary = first(mc_event.mc_trks)
            push!(mc_dset, MCInfo(
                event.header.detector_id,
                event.header.run,
                event.header.frame_index,
                event.header.trigger_counter,
                primary.pos...,
                primary.dir...,
                azimuth(primary.dir),
                zenith(primary.dir),
                primary.t,
                primary.E,
                primary.type,
                primary.type == 81 ? length(mc_event.mc_trks) - 1 : 0,  # 81 is muon bundle
            ))
        end

        for idx in 1:max(msfparams.nfits, length(muons))
            muon = muons[idx]

            push!(reco_dset, FFitResult(
                event.header.detector_id,
                event.header.run,
                event.header.frame_index,
                event.header.trigger_counter,
                muon.pos...,
                muon.dir...,
                azimuth(muon.dir),
                zenith(muon.dir),
                muon.t,
                muon.Q,
                muon.S1,
                muon.S2,
                muon.NDF,
                walltime,
            ))
        end
    end
    close(outfile)

    println("Total number of events processed: $(n.value)")
    println("Number of muon candidates: $(n_muons.value)")
    println("Failed reconstructions: $(n_failed.value)")
end

main()
