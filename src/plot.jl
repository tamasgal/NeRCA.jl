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
