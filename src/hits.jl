struct Hit <: KM3io.AbstractHit
    t::Float64
    tot::Float64
end
Base.isless(lhs::Hit, rhs::Hit) = lhs.t < rhs.t

struct HitL0 <: KM3io.AbstractHit
    channel_id::Int8
    t::Float64
    tot::Float64
    pos::Position{Float64}
    dir::Direction{Float64}
end
Base.isless(lhs::HitL0, rhs::HitL0) = time(lhs) < time(rhs)

abstract type AbstractSpecialHit <: KM3io.AbstractHit end  # TODO: bad naming ;)
abstract type AbstractCoincidenceHit <: AbstractSpecialHit end
abstract type AbstractReducedHit <: AbstractSpecialHit end

struct HitL1 <: AbstractCoincidenceHit
    dom_id::Int32
    hits::Vector{HitL0}
end
struct HitL2 <: AbstractCoincidenceHit
    dom_id::Int32
    hits::Vector{HitL0}
end
HitL1(m::DetectorModule, hits) = HitL1(m.id, hits)
HitL2(m::DetectorModule, hits) = HitL2(m.id, hits)
Base.length(c::AbstractCoincidenceHit) = length(c.hits)
Base.eltype(c::AbstractCoincidenceHit) = HitL0
function Base.show(io::IO, c::AbstractCoincidenceHit)
    times = [time(h) for h in c]
    print(io, "$(c.dom_id) $(minimum(times)) $(length(c))")
end
function Base.iterate(c::AbstractCoincidenceHit, state=1)
    @inbounds state > length(c) ? nothing : (c.hits[state], state+1)
end

struct HitR0 <: AbstractReducedHit
    hit::Hit
    channel_id::Int8
end
Base.isless(lhs::HitR0, rhs::HitR0) = lhs.hit.t < rhs.hit.t

struct HitR1 <: AbstractReducedHit
    dom_id::Int32
    pos::Position{Float64}
    hit::Hit
    n::Int
    weight::Float64
end
Base.isless(lhs::HitR1, rhs::HitR1) = lhs.dom_id == rhs.dom_id ? lhs.hit.t < rhs.hit.t : lhs.dom_id < rhs.dom_id
const HitR2 = HitR1
function HitR1(m::DetectorModule, hits::Vector{HitL0})
    combined_hit = combine(hits)
    h = first(hits)
    count = weight = length(hits)
    HitR1(m.id, h.pos, combined_hit, count, weight)
end

starttime(hit) = time(hit)
endtime(hit) = time(hit) + hit.tot

"""
Combine snapshot and triggered hits to a single hits-vector.

This should be used to transfer the trigger information to the
snapshot hits from a DAQEvent. The triggered hits are a subset
of the snapshot hits.

"""
function combine(snapshot_hits::Vector{KM3io.SnapshotHit}, triggered_hits::Vector{KM3io.TriggeredHit})
    triggermasks = Dict{Tuple{UInt8, Int32, Int32, UInt8}, Int64}()
    for hit ∈ triggered_hits
        triggermasks[(hit.channel_id, hit.dom_id, hit.t, hit.tot)] = hit.trigger_mask
    end
    n = length(snapshot_hits)
    hits = sizehint!(Vector{TriggeredHit}(), n)
    for hit in snapshot_hits
        channel_id = hit.channel_id
        dom_id = hit.dom_id
        t = hit.t
        tot = hit.tot
        triggermask = get(triggermasks, (channel_id, dom_id, t, tot), 0)
        push!(hits, TriggeredHit(dom_id, channel_id, t, tot, triggermask))
    end
    hits
end


"""
Create a `Vector` with hits contributing to `n`-fold coincidences within a time
window of Δt.
"""
function nfoldhits(hits::Vector{T}, Δt, n) where {T<:KM3io.AbstractDAQHit}
    hit_map = modulemap(hits)
    chits = Vector{T}()
    for (dom_id, dom_hits) ∈ hit_map
        bag = Vector{T}()
        push!(bag, dom_hits[1])
        t0 = dom_hits[1].t
        for hit in dom_hits[2:end]
            if hit.t - t0 > Δt
                if length(bag) >= n
                    append!(chits, bag)
                end
                bag = Vector{T}()
            end
            push!(bag, hit)
            t0 = hit.t
        end
    end
    return chits
end


"""
Calculate the multiplicities for a given time window. Two arrays are
are returned, one contains the multiplicities, the second one the IDs
of the coincidence groups.
The hits should be sorted by time and then by dom_id.
"""
function count_multiplicities(hits::Vector{T}, tmax=20) where {T<:KM3io.AbstractHit}
    n = length(hits)
    mtp = ones(Int32, n)
    cid = zeros(Int32, n)
    idx0 = 1
    _mtp = 1
    _cid = idx0
    t0 = hits[idx0].t
    dom_id = hits[idx0].dom_id
    for i in 2:n
        hit = hits[i]
        if hit.dom_id != dom_id
            mtp[idx0:i-1] .= _mtp
            cid[idx0:i-1] .= _cid
            dom_id = hit.dom_id
            t0 = hit.t
            _mtp = 1
            _cid += 1
            idx0 = i
            continue
        end
        Δt = hit.t - t0
        if Δt > tmax
            mtp[idx0:i] .= _mtp
            cid[idx0:i] .= _cid
            _mtp = 0
            _cid += 1
            idx0 = i
            t0 = hit.t
        end
        _mtp += 1
    end
    mtp[idx0:end] .= _mtp
    cid[idx0:end] .= _cid
    mtp, cid
end

"""
Counts the multiplicities and modifies the .multiplicity field of the hits.
Important: the hits have to be sorted by time and then by DOM ID first.
"""
function count_multiplicities!(hits::Vector{KM3io.XCalibratedHit}, tmax=20)
    _mtp = 0
    _cid = 0
    t0 = 0
    dom_id = 0
    hit_buffer = Vector{XCalibratedHit}()

    function process_buffer()
        while !isempty(hit_buffer)
            _hit = pop!(hit_buffer)
            _hit.multiplicity.count = _mtp
            _hit.multiplicity.id = _cid
        end
    end

    function reset()
        _mtp = 1
        _cid += 1
    end

    for hit ∈ hits
        if length(hit_buffer) == 0
            reset()
            push!(hit_buffer, hit)
            t0 = hit.t
            dom_id = hit.dom_id
            continue
        end
        if hit.dom_id != dom_id
            process_buffer()
            push!(hit_buffer, hit)
            t0 = hit.t
            dom_id = hit.dom_id
            reset()
            continue
        end
        Δt = hit.t - t0
        if Δt > tmax
            process_buffer()
            push!(hit_buffer, hit)
            t0 = hit.t
            reset()
        else
            push!(hit_buffer, hit)
            _mtp += 1
        end
    end
    if length(hit_buffer) > 0
        process_buffer()
    end
    return
end


"""
Categorise hits by DU and put them into a dictionary of DU=>Vector{Hit}.

Caveat: this function is not typesafe, only suited for high-level analysis (like plots).
"""
@inline duhits(hits::Vector{T}) where {T<:KM3io.XCalibratedHit} = categorize(:du, hits)


"""
Return a vector of hits with ToT >= `tot`.
"""
function totcut(hits::Vector{T}, tot) where {T<:KM3io.AbstractDAQHit}
    return filter(h->h.tot >= tot, hits)
end


"""
Returns the estimated number of photoelectrons for a given ToT.
"""
function nphes(tot)
    if tot <= 20
        return 0.0
    end
    if tot <= 26
        return 1.0
    end
    if tot < 170
        return 1.0 + (tot - 26)/(1/0.28)
    end
    return 40.0 + (255 - tot)*2.0
end


"""

Creates a map (`Dict{Int32, Vector{T}}`) from a flat `Vector{T}` split up based
on the `dom_id` of each element. A typical use is to split up a vector of hits
by their optical module IDs.

This function is similar to `categorize(:dom_id, Vector{T})` but this method
is completely typesafe.

"""
function modulemap(hits::Vector{T}) where T
    out = Dict{Int32, Vector{T}}()
    for hit ∈ hits
        if !(hit.dom_id ∈ keys(out))
            out[hit.dom_id] = T[]
        end
        push!(out[hit.dom_id], hit)
    end
    out
end

"""
Calibrates hits.
"""
function calibrate(T::Type{HitR1}, det::Detector, hits)
    rhits = sizehint!(Vector{T}(), length(hits))
    for hit ∈ hits
        pmt = det[hit.dom_id][hit.channel_id]
        t = hit.t + pmt.t₀
        push!(rhits, T(hit.dom_id, pmt.pos, Hit(t, hit.tot)))
    end
    rhits
end
function calibrate(T::Type{HitL0}, m::DetectorModule, hits)
    chits = sizehint!(Vector{T}(), length(hits))
    for hit ∈ hits
        pmt = m[hit.channel_id]
        t = hit.t + pmt.t₀
        push!(chits, T(hit.channel_id, t, hit.tot, pmt.pos, pmt.dir))
    end
    chits
end


"""
Combines several hits into a single one by taking the earliest start time,
then latest endtime and a ToT which spans over the time range.
"""
function combine(hits::Vector{T}) where T <: KM3io.AbstractHit
    isempty(hits) && return Hit(0.0, 0.0)

    t1 = starttime(first(hits))
    t2 = endtime(first(hits))

    @inbounds for hit ∈ hits[2:end]
        _t1 = starttime(hit)
        _t2 = endtime(hit)
        if t1 > _t1
            t1 = _t1
        end
        if t2 < _t2
            t2 = _t2
        end
    end
    Hit(t1, t2 - t1)
end


struct L1BuilderParameters
    Δt::Float64
    combine::Bool
end

struct L1Builder
    params::L1BuilderParameters
end

"""
Find coincidences within the time window `Δt` of the initialised `params`. The return
value is a vector of `L1Hit`s.
"""
function (b::L1Builder)(::Type{H}, det::Detector, hits::Vector{T}; combine=false) where {T, H}
    out = H[]
    mm = modulemap(hits)
    for (m, module_hits) ∈ mm
        _findL1!(out, det[m], module_hits, b.params.Δt; combine=combine)
    end
    out
end
(b::L1Builder)(det::Detector, hits) = b(HitL1, det, hits)
function _findL1!(out::Vector{H}, m::DetectorModule, hits, Δt; combine=false) where H <: AbstractSpecialHit
    n = length(hits)
    n < 2 && return out

    chits = sort(calibrate(HitL0, m, hits))

    last_coincidence_idx = 0  # used when combining

    ref_idx = 1  # starting with the first hit, obviously
    idx = 2      # first comparison is the second hit
    while ref_idx <= n
        restart = false
        idx = ref_idx + 1
        while(idx <= n+1 || restart) # n+1 to do the final loop after the last hit
            if idx > n || (time(chits[idx]) - time(chits[ref_idx]) > Δt)
                end_idx = idx - 1
                # check if we have gathered some hits
                if ref_idx != end_idx
                    coincident_hits = [chits[i] for i ∈ ref_idx:end_idx]
                    push!(out, H(m, coincident_hits))
                    ref_idx = end_idx
                end
                restart = true
                break
            end
            idx += 1
        end
        ref_idx += 1
    end
    out
end


struct L2BuilderParameters
    n_hits::Int
    Δt::Float64
    ctmin::Float64
end

struct L2Builder
    params::L2BuilderParameters
end

function (b::L2Builder)(hits)
end


const SLEWS_L1 = SVector(
    +0.00, +0.39, +0.21, -0.59, -1.15,
    -1.59, -1.97, -2.30, -2.56, -2.89,
    -3.12, -3.24, -3.56, -3.69, -4.00,
    -4.10, -4.16, -4.49, -4.71, -4.77,
    -4.81, -4.87, -4.88, -4.83, -5.21,
    -5.06, -5.27, -5.18, -5.24, -5.79,
    -6.78, -6.24
)

function KM3io.slew(h::HitR1)
    h.n > length(SLEWS_L1) && return SLEWS_L1[end]
    SLEWS_L1[h.n]
end
