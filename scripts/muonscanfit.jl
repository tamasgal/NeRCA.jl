#!/usr/bin/env julia
doc = """Muon Scanfit

Usage:
  muonscanfit.jl -a DETECTORFILE -i INPUTFILE -o OUTPUTFILE
  muonscanfit.jl -h | --help
  muonscanfit.jl --version

Options:
  -a DETECTORFILE  The DETX file containing geometry and time calibration.
  -i INPUTFILE     ROOT file containing "online" events.
  -o OUTPUTFILE    The HDF5 file to write the output into.
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
    f = NeRCA.ROOTFile(args["-i"])
    det = KM3io.Detector(args["-a"])

    msfparams = MuonScanfitParameters(;tmaxlocal=18.0, roadwidth=200.0)
    msfit = MuonScanfit(msfparams, det)

    outfile = H5File(args["-o"], "w")
    dset = create_dataset(outfile, "reco/muonscanfit", MuonScanfitResult)

    n_events = 0
    n_muons = 0
    n_failed = 0
    @showprogress 1 for event ∈ f.online.events
        walltime = @elapsed muons = msfit(event.snapshot_hits)

        n_events += 1
        n_muons += length(muons)

        if length(muons) == 0
            n_failed += 1
            continue
        end

        best_muon = first(muons)

        mc_angular_error = NaN
        mc_energy = NaN
        if f.online !== nothing
            neutrino = first(f.offline[event.header.trigger_counter + 1].mc_trks)
            mc_angular_error = rad2deg(angle(best_muon.dir, neutrino.dir))
            mc_energy = neutrino.E
        end

        push!(dset, MuonScanfitResult(
            event.header.detector_id,
            event.header.run,
            event.header.frame_index,
            event.header.trigger_counter,
            best_muon.pos...,
            best_muon.dir...,
            best_muon.t,
            best_muon.Q,
            best_muon.NDF,
            walltime,
            mc_angular_error,
            mc_energy,
        ))
    end
    close(outfile)

    println("Total number of events processed: $n_events")
    println("Number of muon candidates: $n_muons")
    println("Failed reconstructions: $n_failed")
end

main()
