"""
    function single_du_params(track::KM3NeT.Track)

Calculates five parameters to describe a track for a single DU case.
"""
function single_du_params(track::KM3NeT.Track)
    pos = track.pos
    dir = Direction(normalize(track.dir))
    t₀ = track.time
    proj = pos ⋅ dir
    z_closest = (pos.z - dir.z*proj) / (1 - dir.z^2)
    t_closest = t₀ + (z_closest * dir.z - proj)/c
    p_t_closest = pos + c * (t_closest - t₀) * dir
    d_closest = √(p_t_closest[1]^2 + p_t_closest[2]^2)

    d_closest, t_closest, z_closest, dir.z, t₀
end
