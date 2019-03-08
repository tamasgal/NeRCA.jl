"""
    function single_du_params(track::KM3NeT.Track)

Calculates five parameters to describe a track for a single DU case.
"""
function single_du_params(track::KM3NeT.Track)
    c_ns = c / 1e9
    pos = track.pos
    dir = Direction(normalize(track.dir))
    t₀ = track.time
    proj = pos ⋅ dir
    z_closest = (pos.z - dir.z*proj) / (1 - dir.z^2)
    t_closest = t₀ + (z_closest * dir.z - proj)/c_ns
    p_t_closest = pos + c_ns * (t_closest - t₀) * dir
    d_closest = √(p_t_closest[1]^2 + p_t_closest[2]^2)

    d_closest, t_closest, z_closest, dir.z, t₀
end


"""
    function single_du_params(t::KM3NeT.MCTrack)

Calculates five parameters to describe a MC track for a single DU case.
"""
function single_du_params(t::KM3NeT.MCTrack, time)
    single_du_params(KM3NeT.Track([t.dir_x, t.dir_y, t.dir_z], [t.pos_x, t.pos_y, t.pos_z], time))
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


function reco(hits, du, calib)
    chits = KM3NeT.calibrate(hits, calib);
    dhits = filter(h -> h.du == du, chits);
    dthits = filter(h -> h.triggered, dhits);
    sort!(dthits, by = h -> h.t);
    shits = unique(h -> h.dom_id, dthits);
    
    qfunc = KM3NeT.make_quality_function([h.pos.z for h in shits], [h.t for h in shits])
    
    model = Model(with_optimizer(Ipopt.Optimizer))

    register(model, :qfunc, 5, qfunc, autodiff=true)

    @variable(model, -1000 <= d_closest <= 1000, start=100.0)
    @variable(model, t_closest, start=shits[1].t)
    @variable(model, -1000 <= z_closest <= 1000, start=shits[1].pos.z)
    @variable(model, -1 <= dir_z <= 1, start=-0.9)
    @variable(model, t₀, start=shits[1].t)
    
    @NLobjective(model, Min, qfunc(d_closest, t_closest, z_closest, dir_z, t₀))
    
    optimize!(model);
    
    return value(d_closest), value(t_closest), value(z_closest), value(dir_z), value(t₀)
end
