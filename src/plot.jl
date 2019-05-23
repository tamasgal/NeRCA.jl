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
    f(hits::Vector{CalibratedHit}; label="hits", highlight_triggered=false, multiplicities=false, Δt=20, du=0, t₀=0)

Plot recipe to plot simple z-t-plots.
"""
@recipe function f(hits::Vector{CalibratedHit}; label="hits", markersize=4, highlight_triggered=false, multiplicities=false, Δt=20, du=0, t₀=0)
    seriestype := :scatter

    xlabel := "time [ns]"
    ylabel := "z [m]"
    markerstrokewidth := 0

    if du > 0
        hits = filter(h -> h.du == du, hits)
    end

    thits = filter(h -> h.triggered, hits)

    @series begin
        label := label
        if multiplicities
            markersize := count_multiplicities(hits)[1]
            markeralpha := 0.8
        end
        [h.t - t₀ for h in hits], [h.pos.z for h in hits]
    end

    if highlight_triggered
        @series begin
            label := "triggered"
            markersize := 1
            marker := :x
            markerstrokewidth := 1
            [h.t - t₀ for h in thits], [h.pos.z for h in thits]
        end
    end
end
