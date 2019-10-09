"""
    function make_mc_time_converter(event_info::Union{MCEventInfo,DAQEventInfo})

Returns a function which converts MC time to JTE time.
"""
function make_mc_time_converter(event_info::Union{MCEventInfo,DAQEventInfo})
    function time_converter(time)
        return time - (event_info.timestamp * 1e9 + event_info.nanoseconds) + event_info.mc_time
    end
    return time_converter
end


"""
    function cherenkov_origin(pos, t::Track)

Calculate the origin of the Cherenkov photon on a track.
"""
function cherenkov_origin(pos, t::Track)
    θ = acos(1/N_SEAWATER)
    P = project(pos, t) 
    dir = -normalize(t.dir)
    track_distance = pld3(pos, t.pos, t.dir)
    distance = track_distance / tan(θ)
    Position(P + distance*dir)
end
