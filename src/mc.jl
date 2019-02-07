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
function make_cherenkov_calculator(track_pos, track_dir; theta=0.759296, c_water=217445751.79)
    tan_theta = tan(theta)
    sin_theta = sin(theta)
    one_over_c_water = 1 / c_water
    one_over_c = 1 / 2.99792458e8
    function cherenkov_time(pos)
        v = pos - track_pos
        l = dot(v, track_dir)
        p = dot(v, v) - l^2
        if p < 0
            return NaN
        end
        k = sqrt(p)
        a_1 = k / tan_theta
        a_2 = k / sin_theta
        t_c = one_over_c * (l - a_1) + one_over_c_water * a_2
        return t_c * 1e9
    end
    return cherenkov_time
end


"""
    function make_mc_time_converter(event_info::MCEventInfo)

Returns a function which converts MC time to JTE time.
"""
function make_mc_time_converter(event_info::MCEventInfo)
    function time_converter(time)
        return time - (event_info.timestamp * 1e9 + event_info.nanoseconds)  \
               + event_info.mc_time
    end
    return time_converter
end
