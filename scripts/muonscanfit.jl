#!/usr/bin/env julia

using LinearAlgebra
BLAS.set_num_threads(1)

doc = """Muon Scanfit

Usage:
  muonscanfit.jl [options] -a DETECTORFILE -i INPUTFILE -o OUTPUTFILE
  muonscanfit.jl -h | --help
  muonscanfit.jl --version

Options:
  -a DETECTORFILE      The DETX file containing geometry and time calibration.
  -i INPUTFILE         ROOT file containing "online" events.
  -n NUMBER_OF_EVENTS  The number of events to process [default: 999999999].
  -o OUTPUTFILE        The HDF5 file to write the output into.
"""
using DocOpt
const args = docopt(doc)  # fail early ;)

println("Loading libraries...")
using NeRCA
using KM3io
using ProgressMeter


struct MuonScanfitResult
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
    t::Float64
    Q::Float64
    NDF::Int
    walltime::Float64
    mc_angular_error::Float64
    mc_energy::Float64
end


function main()
    f = ROOTFile(args["-i"])
    det = Detector(args["-a"])
    n_events = parse(Int, args["-n"])

    msfparams = MuonScanfitParameters(;tmaxlocal=18.0, roadwidth=200.0, nfits=12)
    msfit = MuonScanfit(msfparams, det)

    outfile = H5File(args["-o"], "w")
    dset = create_dataset(outfile, "reco/muonscanfit", MuonScanfitResult)

    addmeta(dset, msfparams)

    n = 0
    n_muons = 0
    n_failed = 0
    @showprogress 1 for event_idx in 1:min(length(f.online.events), n_events)
        event = f.online.events[event_idx]

        walltime = @elapsed muons = msfit(event.snapshot_hits)

        n += 1
        n_muons += length(muons)

        if length(muons) == 0
            n_failed += 1
            continue
        end

        nu = f.online !== nothing ? first(f.offline[event.header.trigger_counter + 1].mc_trks) : missing

        for idx in 1:max(msfparams.nfits, length(muons))
            muon = muons[idx]
            mc_angular_error = NaN
            mc_energy = NaN
            if !ismissing(nu)
                mc_angular_error = rad2deg(angle(muon.dir, nu.dir))
                mc_energy = nu.E
            end

            push!(dset, MuonScanfitResult(
                event.header.detector_id,
                event.header.run,
                event.header.frame_index,
                event.header.trigger_counter,
                muon.pos...,
                muon.dir...,
                muon.t,
                muon.Q,
                muon.NDF,
                walltime,
                mc_angular_error,
                mc_energy,
            ))
        end
    end
    close(outfile)

    println("Total number of events processed: $n")
    println("Number of muon candidates: $n_muons")
    println("Failed reconstructions: $n_failed")
end

main()
