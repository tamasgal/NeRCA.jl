"""
    function cherenkov_origin(pos, t::KM3io.Track)

Calculate the origin of the Cherenkov photon on a track.
"""
function cherenkov_origin(pos, t::KM3io.Track)
    θ = acos(1/KM3io.INDEX_OF_REFRACTION_WATER)
    P = project(pos, t) 
    dir = -normalize(t.dir)
    track_distance = pld3(pos, t.pos, t.dir)
    distance = track_distance / tan(θ)
    Position(P + distance*dir)
end
