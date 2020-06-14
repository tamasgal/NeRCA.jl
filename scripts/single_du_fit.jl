println("Loading libraries...")
using Distributed

addprocs(4)

@everywhere begin
    using Pkg
    Pkg.activate(".")
    using NeRCA
    using DrWatson
    using ProgressMeter

    struct RecoFile{T}
        filepath::AbstractString
        _fobj

        function RecoFile{T}(filepath::AbstractString) where T
            if isfile(filepath)
                @warn "Reconstruction file '$filepath' present, overwriting..."
            end
            fobj = open(filepath, "w")
            write(fobj, join(["event_id", "du", "Q", "n_triggered_doms", "n_hits", "n_triggered_hits", fieldnames(T)...], ",") * "\n")
            new(filepath, fobj)
        end
    end

    close(f::RecoFile) = close(f._fobj)

    function Base.write(f::RecoFile, event_id, du, Q, n_triggered_doms, n_hits, n_triggered_hits, s::SingleDUParams)
        write(f._fobj, "$event_id,$du,$Q,$n_triggered_doms,$n_hits,$n_triggered_hits,")
        write(f._fobj, join([getfield(s, field) for field in fieldnames(typeof(s))], ","))
        write(f._fobj, "\n")
    end
end

if length(ARGS) < 3
    println("Usage: ./single_du_fit.jl OUTPATH DETX ROOTFILE")
    exit(1)
end


function main()
    println("Starting reconstruction.")

    sparams = SingleDURecoParams(max_iter=200)

    outpath = ARGS[1]
    detx = ARGS[2]
    rootfile = ARGS[3]
    outfile = joinpath(outpath, basename(rootfile) * ".ROyFit.csv")

    mkpath(outpath)
    recofile = RecoFile{SingleDUParams}(outfile)

    calib = Calibration(detx)
    f = NeRCA.OnlineFile(rootfile)

    event_shits = NeRCA.read_snapshot_hits(f)
    event_thits = NeRCA.read_triggered_hits(f)
    n_events = length(event_shits)
    println("$n_events events found")

    @showprogress pmap(zip(1:n_events, event_shits, event_thits)) do (event_id, shits, thits)
    # for (shits, thits) in zip(event_shits, event_thits)
        hits = calibrate(calib, NeRCA.combine(shits, thits))


        triggered_hits = triggered(hits)

        dus = sort(unique(map(h->h.du, hits)))
        triggered_dus = sort(unique(map(h->h.du, triggered_hits)))
        n_dus = length(dus)
        n_triggered_dus = length(triggered_dus)

        for du in dus
            du_hits = filter(h->h.du == du, hits)
            triggered_du_hits = triggered(du_hits)
            n_triggered_hits = length(triggered_du_hits)
            n_triggered_doms = length(unique(h->h.dom_id, triggered_du_hits))
            if n_triggered_hits == 0
                continue
            end
            fit = NeRCA.single_du_fit(du_hits, sparams)
            write(recofile, event_id, du, fit.Q, n_triggered_doms, length(du_hits), n_triggered_hits, fit.sdp)
        end
    end
end


function initialise_recofile(path)
end


main()
