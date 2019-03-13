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
    function make_quality_function(z_positions, times)

Generates a function to be used in an optimiser. The generated function has
the following signature:

    quality_function(d_closest, t_closest, z_closest, dir_z, t₀)
"""
function make_quality_function(hits::Vector{KM3NeT.CalibratedHit})
    z_positions, times = [h.pos.z for h in hits], [h.t for h in hits]
    tots = [h.tot for h in hits]
    charges = map(t -> t <= 26 ? 0 : floor((t - 26) / 7) + 1, tots)
    mean_charge = mean(charges)
    mean_tot = mean(tots)
    print(map(signed, tots))
    function quality_function(d_closest, t_closest, z_closest, dir_z, t₀)
        d_γ, ccalc = make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)
        expected_times = ccalc.(z_positions)
        #= delta_zs = abs.(z_positions .- z_closest) =#
        #= delta_ts = delta_zs ./ maximum(delta_zs) .* 2 =#
        Δts = abs.(times - expected_times)
        return sum(Δts .^2 + tots ./ mean_tot)
    end
    return quality_function
end


function reco(hits::Vector{KM3NeT.CalibratedHit}; print_level=0)
    dthits = filter(h -> h.triggered, hits);
    sort!(dthits, by = h -> h.t);
    shits = unique(h -> h.dom_id, dthits);

    qfunc = KM3NeT.make_quality_function(shits)

    model = Model(with_optimizer(Ipopt.Optimizer, print_level=print_level, tol=1e-3))

    register(model, :qfunc, 5, qfunc, autodiff=true)

    hit_time = shits[1].t

    d_closest_start = 100.0
    t_closest_start = hit_time
    z_closest_start = shits[1].pos.z
    dir_z_start = -0.9
    t₀_start = hit_time

    max_z = maximum([h.pos.z for h in shits]) - 2*38
    min_z = minimum([h.pos.z for h in shits]) + 2*38
    if min_z > max_z
        min_z = abs(max_z - min_z) / 2
    end

    #= d_closest_start = 0.0 =#
    #= t_closest_start = 0 =#
    #= z_closest_start = 10 =#
    #= dir_z_start = 0 =#
    #= t₀_start = 0 =#

    @variable(model, 1 <= d_closest <= 1000, start=d_closest_start)
    @variable(model, hit_time - 1000 <= t_closest <= hit_time + 1000, start=t_closest_start)
    @variable(model, min_z <= z_closest <= max_z, start=z_closest_start)
    @variable(model, -1 <= dir_z <= 1, start=dir_z_start)
    @variable(model, t₀, start=t₀_start)

    @NLobjective(model, Min, qfunc(d_closest, t_closest, z_closest, dir_z, t₀))

    optimize!(model);

    values = (value(d_closest), value(t_closest), value(z_closest), value(dir_z), value(t₀))

    return values, qfunc(values...)/4, model
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
