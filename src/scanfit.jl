struct MuonScanfitParameters
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    directions::Vector{Direction{Float64}}
end
function Base.show(io::IO, mp::MuonScanfit)
   print(io, "$(typeof(mp)) in $(length(mp.directions)) directions")
end

function (mp::MuonScanfit)(hits)
    rhits = calibrate(HitR1, mp.detector, hits)

    for dir ∈ mp.directions
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits)
            rhits[idx] = @set rhit.pos = R * rhit.pos
        end



        # rotate hits back
        for (idx, rhit) ∈ enumerate(rhits)
            rhits[idx] = @set rhit.pos = R \ rhit.pos
        end
    end
    rhits
end

struct Line1Z
    pos::Position{Float64}
    t::Float64
end

struct Hit
    t::Float64
    tot::Float64
end

struct HitR1
    dom_id::Int32
    pos::Position{Float64}
    hit::Hit
end
Base.isless(lhs::HitR1, rhs::HitR1) = lhs.dom_id == rhs.dom_id ? lhs.t < rhs.t : lhs.dom_id < rhs.dom_id

starttime(hit) = hit.t
endtime(hit) = hit.t + hit.tot

function calibrate(T::Type{HitR1}, det::Detector, hits)
    rhits = sizehint!(Vector{T}(), length(hits))
    for hit ∈ hits
        pos = det[hit.dom_id][hit.channel_id].pos
        push!(rhits, T(hit.dom_id, pos, Hit(hit.t, hit.tot)))
    end
    rhits
end

"""
Combines several hits into a single one by taking the earlierst start time,
then latest endtime and a ToT which spans over the time range.
"""
function combine(hits::Vector{Hit})
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

struct L2BuilderParameters
    n_hits::Int
    Δt::Float64
    ctmin::Float64
end

struct L2Builder
    params::L2BuilderParameters
end

"""

Creates a map (`Dict{Int32, Vector{T}}`) from a flat `Vector{T}` split up
based on the `dom_id`. A typical use is to split up a vector of hits
by the optical module ID.

This function is similar to `categorize(:dom_id, Vector{T})` but typesafe.

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
