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
    function svdfit(hits::Vector{CalibratedHit})

Uses SVD to do a fast and dirty track prefit. Provide hits with a multiplicity
of at least 2.
"""
function svdfit(hits::Vector{CalibratedHit})
    t₀ = hits[div(length(hits), 2)].t
    pos, dir = svdfit(matrix([h.pos for h in hits]))
    if (last(hits).pos - first(hits).pos) ⋅ dir  < 0.0
        dir *= -1
    end
    return Track(dir, pos, t₀)
end


"""
    function prefit(hits::Vector{CalibratedHit})

Performs the prefit algorithm which was used in DUMAND II.
"""
function prefit(hits::Vector{CalibratedHit})
    N = length(hits)
    D = 0.0
    dir = [0.0, 0.0, 0.0]
    pos = [0.0, 0.0, 0.0]
    pes = [max(1, (h.tot - 26.3) / 4.5) for h in hits]
    # pes = [h.tot for h in hits]
    # pes = [h.multiplicity.count for h in hits]
    @inbounds for i ∈ 1:N
        @inbounds for k ∈ 1:N
            if i == k
                continue
            end
            t_i = hits[i].t
            t_k = hits[k].t
            q_ik = pes[i] * pes[k]
            t_ik = t_i * t_k
            pos += q_ik * (hits[i].pos*(t_k^2 - t_ik) + hits[k].pos*(t_i^2 - t_ik))
            dir += q_ik * (hits[i].pos - hits[k].pos) * (t_i - t_k)
            D += q_ik * (t_i - t_k)^2
        end
    end
    # dir = normalize(dir/D)
    return NeRCA.Track(dir/D, pos/D, 0)
end


"""
    function make_cherenkov_calculator(track::Track; v=2.99792458e8)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(track::Track; v=2.99792458e8, n=N_SEAWATER)
    c_medium = c_0/n
    β = v/c_0
    θ = acos(min(1/(β*n), 1))
    θ′ = π - θ
    track_dir = normalize(track.dir)
    t₀ = track.time

    pos::Position -> begin

        # minimal distance between hit and track
        track_distance = pld3(track.pos, pos, track_dir)
        # the path the Cherenkov photon has to travel
        cherenkov_path_length = track_distance / sin(θ)
        distance_vector = pos - track.pos
        # distance between particle position and hit
        distance = norm(distance_vector)

        t_cherenkov = cherenkov_path_length / c_medium*1e9

        η = angle_between(track_dir, distance_vector)
        χ = θ - η
        particle_travel_path_length = distance / sin(θ′) * sin(χ)
        t = t₀ + particle_travel_path_length / v*1e9 + t_cherenkov
        return t
    end
end


"""
    function make_cherenkov_calculator(track::Track, event_info::Union{MCEventInfo,DAQEventInfo}; v=2.99792458e8)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(track::Track, event_info::Union{MCEventInfo,DAQEventInfo}; v=2.99792458e8)
    jte_time = make_mc_time_converter(event_info)(track.time)
    make_cherenkov_calculator(Track(track.dir, track.pos, jte_time), v=v)
end


"""
    function make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀; n=N_SEAWATER)
    c_ns = c_0 / 1e9
    d_γ(z) = n/√(n^2 - 1) * √(d_closest^2 + (z-z_closest)^2 * (1 - dir_z^2))
    t(z) = (t₀) + 1/c_ns * ((z - z_closest)*dir_z + (n^2 - 1)/n * d_γ(z))
    d_γ, t
end

"""
    function make_cherenkov_calculator(sdp::SingleDUParams)

Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(sdp::SingleDUParams; n=N_SEAWATER)
    c_ns = c_0 / 1e9
    d_γ(z) = n/√(n^2 - 1) * √(sdp.d^2 + (z-sdp.z)^2 * (1 - sdp.dz^2))
    t(z) = (sdp.t₀) + 1/c_ns * ((z - sdp.z)*sdp.dz + (n^2 - 1)/n * d_γ(z))
    d_γ, t
end



"""
    function single_du_params(track::NeRCA.Track)

Calculates five parameters to describe a track for a single DU case.
"""
function single_du_params(track::NeRCA.Track)
    c_ns = c_0 / 1e9
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
    function single_du_params(t::NeRCA.MCTrack)

Calculates five parameters to describe a MC track for a single DU case.
"""
function single_du_params(t::NeRCA.MCTrack, time)
    single_du_params(NeRCA.Track([t.dir_x, t.dir_y, t.dir_z], [t.pos_x, t.pos_y, t.pos_z], time))
end


# TODO: delete this function
function cherenkov_plausible(Δt, Δz; time_extra=10, n=N_SEAWATER)
    Δt < Δz * n / NeRCA.c_0*1e9 + time_extra
end


function create_hit_pool(hits::T) where T<:Vector{NeRCA.CalibratedHit}
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


function select_hits(hits, hit_pool; Δt=10, Δz=9, new_hits=nothing)
    intervals = Dict{Int16, Interval{Float64}}()
    hits₀ = NeRCA.create_hit_pool(hits)
    
    extended = Vector{CalibratedHit}()
    
    function time_interval(t₀, tₜ, Δfloor)
        t₋ = tₜ - Δfloor*Δt
        t₊ = max(t₀ + Δfloor*Δz*NeRCA.N_SEAWATER/NeRCA.c_0*1e9 + Δt, tₜ + Δfloor*Δt)
        @interval(t₋, t₊)
    end
    
    function add_if_in_interval(floor, interval)
        # println("    searching for hit candidate on floor $floor")
        if haskey(hit_pool, floor)
            floor_hits = hit_pool[floor]
            # println("       $(length(floor_hits)) hits to check")
            for (iₕ, hit) ∈ enumerate(floor_hits)
                # println("           -> hit time: $(hit.t)")
                if hit.t ∈ interval
                    # println("           --> in interval!")
                    if  haskey(hits₀, floor) && hit.t > hits₀[floor][1].t
                        # println("             !! got a hit already which is earlier")
                        # println("             -- removing consecutive hits")
                        deleteat!(floor_hits, iₕ:length(floor_hits))
                        return
                    end
                    # println("             ++ adding hit")
                    push!(extended, hit)
                    # println("             -- removing consecutive hits")
                    deleteat!(floor_hits, iₕ:length(floor_hits))
                    return
                end
            end
        end
    end
    
    hit_seeds = new_hits == nothing ? hits : new_hits
    
    for hit in hit_seeds
        t₀ = hit.t
        floor = hit.floor
        # println("Checking floor $floor")
        for (floor₋, floor₋₋, floor₊, floor₊₊) in [(floor-1, floor-2, floor+1, floor+2), (floor+1, floor+2, floor-1, floor-2)]
            if haskey(hits₀, floor₋)
                t₋ = first(hits₀[floor₋]).t
                if haskey(hits₀, floor₋₋)
                    # println(" => Got hit on $(floor₋₋) and $(floor₋)")
                    t₋₋ = first(hits₀[floor₋₋]).t
                    t₊ = t₋ + 2(t₋ - t₋₋)
                    t₊₊ = t₋ + 3(t₋ - t₋₋)
                    interval₊ = time_interval(t₀, t₊, 1)
                    interval₊₊ = time_interval(t₀, t₊₊, 2)
                    # println("    interval+: $interval₊")
                    # println("    interval++: $interval₊₊")
                    #interval₊ = @interval(minimum([t₀, t₋, t₋₋]) - Δt, maximum([t₀, t₋, t₋₋]) + Δt)
                    #interval₊₊ = @interval(minimum([t₀, t₋, t₋₋]) - 2Δt, maximum([t₀, t₋, t₋₋]) + 2Δt)
                    add_if_in_interval(floor₊, interval₊)
                    add_if_in_interval(floor₊₊, interval₊₊)
                else
                    t₊ = t₀ + (t₀ - t₋)
                    t₊₊ = t₀ + 2(t₀ - t₋)
                    interval₊ = time_interval(t₀, t₊, 1)
                    interval₊₊ = time_interval(t₀, t₊₊, 2)
                    add_if_in_interval(floor₊, interval₊)
                    add_if_in_interval(floor₊₊, interval₊₊)
                end
            elseif haskey(hits₀, floor₋₋)
                t₋₋ = first(hits₀[floor₋₋]).t
                t₊ = t₀ + (t₀ - t₋₋)/2
                t₊₊ = t₀ + 2(t₀ - t₋₋)/2
                interval₊ = time_interval(t₀, t₊, 1)
                interval₊₊ = time_interval(t₀, t₊₊, 2)
                add_if_in_interval(floor₊, interval₊)
                add_if_in_interval(floor₊₊, interval₊₊)
            end
        end
        
    end
    
    if length(extended) == 0
        return hits
    end
    extension = unique(h->h.floor, sort(vcat(hits, extended), by=h->h.t))
    return select_hits(extension, hit_pool; Δt=Δt, Δz=Δz, new_hits=extended)
end


struct SingleDUMinimiser <: Function
    hit_positions::Vector{Position}
    z_positions::Vector{Float64}
    times::Vector{Float64}
    pmt_directions::Vector{Direction}
    multiplicities::Vector{Int}
    max_multiplicity::Float64
    nphes::Vector{Float64}
    average_coinc_tots::Vector{Float64}
    x::Float64
    y::Float64
end

function SingleDUMinimiser(hits::Vector{CalibratedHit}, triggered_hits::Vector{CalibratedHit})
    n = length(hits)
    n_triggered = length(triggered_hits)
    hit_positions = Vector{Position}()
    x_positions = Vector{Float64}()
    y_positions = Vector{Float64}()
    z_positions = Vector{Float64}()
    times = Vector{Float64}()
    multiplicities = Vector{Int32}()
    pmt_directions = Vector{Direction}()
    sizehint!(z_positions, n)
    sizehint!(times, n)
    sizehint!(multiplicities, n)
    sizehint!(pmt_directions, n_triggered)
    @inbounds for i ∈ 1:n
        hit = hits[i]
        push!(hit_positions, hit.pos)
        push!(x_positions, hit.pos.x)
        push!(y_positions, hit.pos.y)
        push!(z_positions, hit.pos.z)
        push!(times, hit.t)
        push!(multiplicities, hit.multiplicity.count)
        push!(pmt_directions, hit.dir)
    end
    # for i ∈ 1:n_triggered
    #     push!(pmt_directions, triggered_hits[i].dir)
    # end
    max_multiplicity = maximum(multiplicities)
    nphe = Vector{Float32}()
    average_coinc_tots = Vector{Float32}()
    for hit in hits
        coinc_hits = filter(h->h.multiplicity.id == hit.multiplicity.id, hits)
        average_coinc_tot = mean([h.tot for h in coinc_hits])
        push!(average_coinc_tots, average_coinc_tot)
        estimated_nphes = sum([nphes(h.tot) for h in coinc_hits])
        push!(nphe, estimated_nphes)
    end
    SingleDUMinimiser(hit_positions, z_positions, times, pmt_directions, multiplicities, max_multiplicity, nphe, average_coinc_tots, mean(x_positions), mean(y_positions))
end


"""
    function (s::SingleDUMinimiser)(d_closest, t_closest, z_closest, dir_z, ϕ, t₀)

The quality function to be minimised when performing the single DU fit.
`d_closest` is the closest distance between the track (starting at `t₀` and
the DU at time `t_closest` with the z-coordinate `z_closest`. The direction
vector is `(0, 0, dir_z)` and `ϕ` the azimuth angle.

"""
function (s::SingleDUMinimiser)(d_closest, t_closest, z_closest, dir_z, ϕ, t₀)
    n = length(s.times)

    d_γ, ccalc = make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀)
    expected_times = ccalc.(s.z_positions)
    # _pos = cartesian(ϕ, π/2; r=d_closest)  # position of the track relative to (0,0,0)
    # track_pos = Position(_pos.x + s.x, _pos.y + s.y, z_closest)
    #
    # track_dir = cartesian(ϕ, acos(dir_z))
    # track = Track(track_dir, track_pos, t₀)

    Q = 0.0
    @inbounds for i ∈ 1:n
        t = s.times[i]
        z = s.z_positions[i]
        m = s.multiplicities[i]
        t_exp = ccalc(z)
        Δt = abs(t - t_exp)
        zenith_acceptance = 1 - NeRCA.zenith(Direction(0,0,1))/π
        photon_distance = d_γ(z)

        # pos_cherenkov = cherenkov_origin(s.hit_positions[i], track)
        # photon_dir = s.hit_positions[i] - pos_cherenkov
        # ξ = angle_between(s.pmt_directions[i], -photon_dir)

        pmtᵩ = azimuth(-s.pmt_directions[i])
        # Δϕ = abs(pmtᵩ - ϕ) #% 2π
        # Δϕ = abs(pmtᵩ - ϕ)
        ξ = (ϕ - pmtᵩ - deg2rad(42))^2
        # Δz = z - z_closest 
        # if dir_z < 0  # downgoing
        #     
        # else
        # end

        
        Q += Δt^2 + photon_distance * s.nphes[i] / 72.0 / 10 * zenith_acceptance + ξ * 10 #+ ξ * 10
    end

    Δts = abs.(s.times - expected_times)
    return Q
end

struct MultiDUMinimiser <: Function
    hits::Vector{CalibratedHit}
end

function (m::MultiDUMinimiser)(x, y, z, θ, ϕ, t₀)
    pos = Position(x, y, z)
    # dir = Direction(cos(θ)*cos(ϕ), cos(θ)*sin(ϕ), sin(θ))
    dir = Direction(ϕ, θ)
    track = Track(dir, pos, t₀)

    # ccalc = make_cherenkov_calculator(track, v=v*1e9)
    ccalc = make_cherenkov_calculator(track)

    Q = 0.0
    for hit in m.hits
        t_exp = ccalc(hit.pos)
        Δt = abs(hit.t - t_exp)
        # println("expected: $(t_exp), got: $(hit.t), delta: $(Δt)")
        Q += (Δt - 2)^2
    end
    Q
end

function multi_du_fit(prefit, hits; print_level=0)
    ϕ_start = NeRCA.azimuth(prefit.dir)
    θ_start = NeRCA.zenith(prefit.dir)
    t₀_start = prefit.time
    pos = prefit.pos

    m = NeRCA.MultiDUMinimiser(hits)

    model = Model(
        with_optimizer(Ipopt.Optimizer,
                       print_level=print_level,
                       max_iter=200,
                       tol=1e-3)
       )
    register(model, :qfunc, 6, m, autodiff=true)

    Δpos = 50

    @variable(model, pos.x - Δpos <= x <= pos.x + Δpos, start=pos.x)
    @variable(model, pos.y - Δpos <= y <= pos.y + Δpos, start=pos.y)
    @variable(model, pos.z - Δpos <= z <= pos.z + Δpos, start=pos.z)
    @variable(model, -1.0 <= θ <= 1.0, start=θ_start)
    @variable(model, -3*2π <= ϕ <= 3*2π, start=ϕ_start)
    @variable(model, t₀, start=t₀_start)

    @NLobjective(model, Min, qfunc(x, y, z, θ, ϕ, t₀))

    optimize!(model);

    return Track(Direction(value(ϕ), value(θ)), Position(value(x), value(y), value(z)), value(t₀))
end


@with_kw struct SingleDURecoParams
    floor_distance::Int = 9
    Δt_extra::Int = 10
    data_type::AbstractString = "mupage"
    discard_worst_hit::Bool = false
end

DrWatson.default_prefix(s::SingleDURecoParams) = "SingleDUReco"


function single_du_fit(du_hits::Vector{NeRCA.CalibratedHit}, par::SingleDURecoParams; print_level=0)
    sort!(du_hits, by = h -> h.t)
    sort!(du_hits, by=h->h.dom_id)
    count_multiplicities!(du_hits)

    hit_pool = create_hit_pool(du_hits)

    hits₀ = unique(h->h.floor, nfoldhits(du_hits, 10, 2))
    if length(hits₀) < 3
        hits₀ = unique(h->h.floor, triggered(du_hits))
    end
    shits = select_hits(hits₀, hit_pool; Δt = par.Δt_extra, Δz = par.floor_distance)

    qfunc = SingleDUMinimiser(shits, filter(h->h.triggered, du_hits))

    model = Model(with_optimizer(Ipopt.Optimizer, print_level=print_level, tol=1e-3))

    register(model, :qfunc, 6, qfunc, autodiff=true)

    brightest_floor = NeRCA.most_frequent(h -> h.floor, du_hits)

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

    d_closest_start = 10.0
    t_closest_start = hit_time
    dir_z_start = -0.9
    ϕ_start = π
    t₀_start = hit_time

    max_z = maximum([h.pos.z for h in shits]) - 2*13
    min_z = minimum([h.pos.z for h in shits]) + 2*13
    if min_z > max_z
        min_z = abs(max_z - min_z) / 2
    end

    @variable(model, 1 <= d_closest <= 100, start=d_closest_start)
    @variable(model, hit_time - 1000 <= t_closest <= hit_time + 1000, start=t_closest_start)
    @variable(model, min_z <= z_closest <= max_z, start=z_closest_start)
    @variable(model, -0.999 <= dir_z <= .999, start=dir_z_start)
    @variable(model, -2π <= ϕ <= 4π, start=ϕ_start)
    @variable(model, t₀, start=t₀_start)

    @NLobjective(model, Min, qfunc(d_closest, t_closest, z_closest, dir_z, ϕ, t₀))

    optimize!(model);
    values = map(value, [d_closest, t_closest, z_closest, dir_z, ϕ, t₀])
    Q₀ = qfunc(values...)/length(shits)

    d_γ, ccalc = make_cherenkov_calculator(map(value, [d_closest, t_closest, z_closest, dir_z, t₀])...) 
    if par.discard_worst_hit && length(shits) > 3
        Δts = [abs(ccalc(h.pos.z) - h.t) for h in shits]
        idx = argmax(Δts)
        deleteat!(Δts, idx)
        deleteat!(shits, idx)

        qfunc = SingleDUMinimiser(shits, filter(h->h.triggered, du_hits))
        qfunc_sym = gensym()
        values_try = map(value, [d_closest, t_closest, z_closest, dir_z, ϕ, t₀])
        register(model, qfunc_sym, 6, qfunc, autodiff=true)
        set_NL_objective(model, MOI.MIN_SENSE, :($qfunc_sym($(d_closest), $(t_closest), $(z_closest), $(dir_z), $(ϕ), $(t₀))))
        optimize!(model)

        Q₀ = qfunc(values...)/length(shits)
    end

    values = map(value, [d_closest, t_closest, z_closest, dir_z, ϕ, t₀])
    sdp = SingleDUParams(values...)
    sdp_initial = SingleDUParams(d_closest_start, t_closest_start, z_closest_start, dir_z_start, ϕ_start, t₀_start)

    return ROyFit(sdp, sdp_initial, qfunc(values...)/length(shits), shits, model)
end
