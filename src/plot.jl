"""
    f(dom_positions, track::KM3NeT.Track)

Plot recipe for z-t-plots.
"""
@recipe function f(dom_positions, track::KM3NeT.Track)
    # layout := @layout [a b]
    layout := (2, 2)

    seriestype := :scatter

    xlabel := "time [ns]"
    ylabel := "z [m]"

    @series begin
        subplot := 1
        ccalc = KM3NeT.make_cherenkov_calculator(track)
        ccalc.(dom_positions), [p.z for p in dom_positions]
    end

    xlabel := "x [m]"
    ylabel := "y [m]"

    @series begin
        subplot := 2
        label := "DOMs"
        [p.x for p in dom_positions], [p.y for p in dom_positions]
        [track.pos.x], [track.pos.y]
    end

    @series begin
        subplot := 2
        label := "Track"
        [track.pos.x], [track.pos.y]
    end

    @series begin
        subplot := 3
    end
end


"""
    f(hits::Vector{CalibratedHit}; triggered_only=false, du=0)

Plot recipe to plot simple z-t-plots.
"""
@recipe function f(hits::Vector{CalibratedHit}; triggered_only=false, du=0)
    seriestype := :scatter

    xlabel := "time [ns]"
    ylabel := "z [m]"
    markerstrokewidth := 0

    if du > 0
        hits = filter(h -> h.du == du, hits)
    end

    thits = filter(h -> h.triggered, hits)
    t₀ = minimum([h.t for h in thits])

    if !triggered_only
        @series begin
            label := "hits"
            [h.t - t₀ for h in hits], [h.pos.z for h in hits]
        end
    end

    @series begin
        label := "triggered hits"
        [h.t - t₀ for h in thits], [h.pos.z for h in thits]
    end

end
