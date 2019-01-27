"""
    function mc_run_id(fname::AbstractString)

Generate a (unique) run ID for a given filename.
"""
function mc_run_id(fname::AbstractString)
    bname = basename(fname)
    if contains(bname, "_muatm")
        s = split(split(bname, "_muatm")[2], ".")[1]
        energy_cut, run = split(s, "T")
        return parse(Int, energy_cut) * 10000 + parse(Int, run)
    end
    if contains(bname, "_numuCC_")
        run = split(split(bname, "_numuCC_")[2], ".")[1]
        return parse(Int, run)
    end
    error("Don't know how to generate a proper run ID for '$bname'.")
end


"""
    function make_cherenkov_calculator(track_pos, track_dir)

Returns a function which calculates the arrival time of a Cherenkov photon
emitted at a given position.
"""
function make_cherenkov_calculator(track_pos, track_dir, theta=0.759296, c_water=2.174458)
    function cherenkov_time(pos)
        v = pos - track_pos
        l = dot(v, track_dir)
        k = dot(v, v) - l^2
        a_1 = k / tan(theta)
        a_2 = k / sin(theta)
        t_c = 1 / 2.99792458e8 * (l - a_1) + 1 / c_water * a_2
        return t_c * 1e9
    end
    return cherenkov_time
end
