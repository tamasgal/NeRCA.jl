struct SingleDUParams
    d::Float64
    t::Float64
    z::Float64
    dz::Float64
    ϕ::Float64
    t₀::Float64
end

struct ROyFit
    sdp::SingleDUParams
    sdp_initial::SingleDUParams
    Q::Float64
    selected_hits::Vector{CalibratedHit}
    model::Model
end


"""
    function make_cherenkov_calculator(track::Track; theta=0.759296, c_water=217445751.79)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(track::Track; theta=0.759296, c_water=217445751.79)
    tan_theta = tan(theta)
    sin_theta = sin(theta)
    one_over_c_water = 1 / c_water
    one_over_c = 1 / 2.99792458e8
    pos::Position -> begin
        v = pos - track.pos
        l = dot(v, normalize(track.dir))
        p = dot(v, v) - l^2
        if p < 0
            return NaN
        end
        k = sqrt(p)
        a_1 = k / tan_theta
        a_2 = k / sin_theta
        t_c = one_over_c * (l - a_1) + one_over_c_water * a_2
        return t_c * 1e9 + track.time
    end
end


"""
    function make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)
    c_ns = c / 1e9
    d_γ(z) = n_water/√(n_water^2 - 1) * √(d_closest^2 + (z-z_closest)^2 * (1 - dir_z^2))
    t(z) = (t₀) + 1/c_ns * ((z - z_closest)*dir_z + (n_water^2 - 1)/n_water * d_γ(z))
    d_γ, t
end

"""
    function make_cherenkov_calculator(sdp::SingleDUParams)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(sdp::SingleDUParams)
    c_ns = c / 1e9
    d_γ(z) = n_water/√(n_water^2 - 1) * √(sdp.d^2 + (z-sdp.z)^2 * (1 - sdp.dz^2))
    t(z) = (sdp.t₀) + 1/c_ns * ((z - sdp.z)*sdp.dz + (n_water^2 - 1)/n_water * d_γ(z))
    d_γ, t
end



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


function cherenkov_plausible(Δt, Δz, time_extra=10)
    Δt < Δz * KM3NeT.n_water / KM3NeT.c*1e9 + time_extra
end


function expand_hits!(hits::T, hit_pool::Dict{Int, T}; max_floor_distance=2) where T<:Vector{KM3NeT.CalibratedHit}
    n = length(hits)
    expanded = false
    for i in 1:n
        hit = hits[i]
        time = hit.t
        floor = hit.floor
        z = hit.pos.z

        for Δfloor ∈ -max_floor_distance:max_floor_distance
            _hits = get(hit_pool, floor + Δfloor, T())
            for _hit ∈ _hits
                Δz = abs(z - _hit.pos.z)
                Δt = abs(time - _hit.t)
                if cherenkov_plausible(Δt, Δz)
                    push!(hits, _hit)
                    expanded = true
                end
            end
        end
    end
end


function create_hit_pool(hits::T) where T<:Vector{KM3NeT.CalibratedHit}
    hit_pool = Dict{Int, T}()
    for hit in hits
        if !haskey(hit_pool, hit.floor)
            hit_pool[hit.floor] = T()
            sizehint!(hit_pool[hit.floor], length(hits))
        end
        push!(hit_pool[hit.floor], hit)
    end
    hit_pool
end


function select_hits(hits::T, hit_pool::Dict{Int, T}) where T<:Vector{KM3NeT.CalibratedHit}

    thits = filter(h -> h.triggered, hits)
    shits = unique(h -> h.dom_id, thits)

    expand_hits!(shits, hit_pool)
    expand_hits!(shits, hit_pool)

    shits = unique(h -> h.dom_id, shits)

    return shits
end


struct SingleDUMinimiser <: Function
    z_positions::Vector{Float64}
    times::Vector{Float64}
    pmt_directions::Vector{Direction}
end


function SingleDUMinimiser(hits::Vector{CalibratedHit})
    n = length(hits)
    z_positions = Vector{Float64}()
    times = Vector{Float64}()
    pmt_directions = Vector{Direction}()
    sizehint!(z_positions, n)
    sizehint!(times, n)
    sizehint!(pmt_directions, n)
    for i ∈ 1:n
        hit = hits[i]
        push!(z_positions, hit.pos.z)
        push!(times, hit.t)
        push!(pmt_directions, hit.dir)
    end
    SingleDUMinimiser(z_positions, times, pmt_directions)
end


function (s::SingleDUMinimiser)(d_closest, t_closest, z_closest, dir_z, ϕ, t₀)
    d_γ, ccalc = make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)
    expected_times = ccalc.(s.z_positions)
    Δts = abs.(s.times - expected_times) 
    Δϕs = filter(!isnan, azimuth.(s.pmt_directions)) .- ϕ
    return sum(Δts .^2) + sum(Δϕs.^2)/length(Δϕs)
    # return sum(Δts .^2)
end


function reco(du_hits::Vector{KM3NeT.CalibratedHit}; print_level=0)
    sort!(du_hits, by = h -> h.t)

    hit_pool = create_hit_pool(du_hits)
    shits = select_hits(du_hits, hit_pool)

    qfunc = SingleDUMinimiser(shits)

    model = Model(with_optimizer(Ipopt.Optimizer, print_level=print_level, tol=1e-3))

    register(model, :qfunc, 6, qfunc, autodiff=true)

    brightest_floor = 0
    floor_hits = 0
    for floor in keys(hit_pool)
        _floor_hits = length(hit_pool[floor])
        if _floor_hits > floor_hits
            floor_hits = _floor_hits
            brightest_floor = floor
        end
    end
    println("Brightest floor: $(brightest_floor)")
    hits_on_brightest_floor = filter(h -> h.floor == brightest_floor, du_hits)
    thits_on_brightest_floor = filter(h -> h.triggered, hits_on_brightest_floor)
    if length(thits_on_brightest_floor) == 0
        z_closest_start = hits_on_brightest_floor[1].pos.z
        hit_time = mean([h.t for h in hits_on_brightest_floor])
    else
        closest_hit = thits_on_brightest_floor[1]
        z_closest_start = closest_hit.pos.z
        hit_time = closest_hit.t
    end

    d_closest_start = 50.0
    t_closest_start = hit_time
    dir_z_start = -0.9
    ϕ_start = π
    t₀_start = hit_time

    max_z = maximum([h.pos.z for h in shits]) - 2*38
    min_z = minimum([h.pos.z for h in shits]) + 2*38
    if min_z > max_z
        min_z = abs(max_z - min_z) / 2
    end

    @variable(model, 1 <= d_closest <= 1000, start=d_closest_start)
    @variable(model, hit_time - 1000 <= t_closest <= hit_time + 1000, start=t_closest_start)
    @variable(model, min_z <= z_closest <= max_z, start=z_closest_start)
    @variable(model, -1 <= dir_z <= 1, start=dir_z_start)
    @variable(model, 0 <= ϕ <= 2π, start=ϕ_start)
    @variable(model, t₀, start=t₀_start)

    @NLobjective(model, Min, qfunc(d_closest, t_closest, z_closest, dir_z, ϕ, t₀))

    optimize!(model);

    values = value(d_closest), value(t_closest), value(z_closest), value(dir_z), value(ϕ), value(t₀)
    sdp = SingleDUParams(values...)
    sdp_initial = SingleDUParams(d_closest_start, t_closest_start, z_closest_start, dir_z_start, ϕ_start, t₀_start)

    return ROyFit(sdp, sdp_initial, qfunc(values...)/length(shits), shits, model)
end
