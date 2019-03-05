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


"""
    function make_quality_function(positions, times)

Generates a function to be used in an optimiser. The generated function has
the following signature:

    quality_function(d_closest, t_closest, z_closest, dir_z, t₀)
"""
function make_quality_function(positions, times)
    function quality_function(d_closest, t_closest, z_closest, dir_z, t₀)
        ccalc = make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)
        expected_times = ccalc.(positions)
        return sum((times - expected_times).^2)
    end
    return quality_function
end
