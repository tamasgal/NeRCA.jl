#!/usr/bin/env julia
using KM3NeT
using Plots
using PlotThemes
using Dates
using Measures
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

        hits = calibrate(event.hits, calib)
        triggered_hits = filter(h->h.triggered, hits)
        n_dus = length(unique(h->h.du, triggered_hits))
        n_doms = length(unique(h->h.dom_id, triggered_hits))

        if n_doms < 4
            continue
        end

        triggered_hits = filter(h->h.triggered, hits)
        dus = sort(unique(map(h->h.du, hits)))
        colours = palette(:default)
        plot()
        Q = []
        for (idx, du) in enumerate(dus)
            du_hits = filter(h->h.du == du, hits)
            if length(triggered(du_hits))== 0
                println("No triggered hits")
                continue
            end
            fit = KM3NeT.single_du_fit(du_hits)
            push!(Q, fit.Q)
            plot!(du_hits, fit, markercolor=colours[idx], label="DU $(du)", max_z=calib.max_z)
        end
        if sum(Q) < 200 && n_doms > 12 && n_dus > 1
            println("Plotting...")
            fit_params = "ROy live reconstruction (combined single line): Q=$([round(_Q,digits=2) for _Q in Q])"
            event_params = "Det ID $(event.det_id), Run $(event.run_id), FrameIndex $(event.timeslice_id), TriggerCounter $(event.trigger_counter), Overlays $(event.overlays)"
            time_params = "$(unix2datetime(event.timestamp)) UTC"
            trigger_params = "Trigger: $(is_mxshower(event) ? "MX " : "")$(is_3dmuon(event) ? "3DM " : "")$(is_3dshower(event) ? "3DS " : "")"
            time_params = "$(unix2datetime(event.timestamp)) UTC"

            plot!(title="$(fit_params)\n$(event_params), $(trigger_params)\n$(time_params)", titlefontsize=8, margin=5mm)

            savefig("ztplot.png")
        end
    end
end


main()
