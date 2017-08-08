__precompile__()

module KM3NeT

using DataStructures
using StaticArrays
using HDF5

import Base: +, -, *

export
    Position, Direction,
    EventInfo, Track, Hit, CalibratedHit, RawHit, McHit, TimesliceHit,
    calibrate, triggered, svdfit, nfoldhits,
    read_indices, read_hits, read_tracks, read_calibration, read_event_info,
    svdfit, matrix, rows


include("types.jl")
include("io.jl")
include("rba.jl")


"""
    function calibrate(hits::Vector{RawHit}, calibration::Calibration)

Apply geometry and time calibration to given hits.
"""
function calibrate(hits::Vector{RawHit}, calibration::Calibration)
    calibrated_hits = Vector{CalibratedHit}()
    for hit in hits
        dom_id = hit.dom_id
        t = hit.t
        channel_id = hit.channel_id
        triggered = hit.triggered
        tot = hit.tot
        pos = calibration.pos[dom_id][channel_id+1]
        dir = calibration.dir[dom_id][channel_id+1]
        t0 = calibration.t0[dom_id][channel_id+1]
        du = calibration.du[dom_id]
        floor = calibration.floor[dom_id]
        c_hit = CalibratedHit(channel_id, dom_id, du, floor, t, tot,
                              triggered, pos, dir, t0)
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


# Hits
"""
    function triggered(hits::Vector{T}) where {T<:Hit}

Return a `Vector` of triggered hits.
"""
triggered(hits::Vector{T}) where {T<:Hit} = filter(h->h.triggered, hits)


"""
    function nfoldhits(hits::Vector{T}, Δt, n) where {T<:Hit}

Create a `Vector` with hits contributing to `n`-fold coincidences within a time
window of Δt.
"""
function nfoldhits(hits::Vector{T}, Δt, n) where {T<:Hit}
    hit_map = DefaultDict{Integer}{Vector{T}}(() -> T[])
    for hit ∈ sort(hits)
        push!(hit_map[hit.dom_id], hit)
    end
    chits = Vector{T}()
    for (dom_id, dom_hits) ∈ hit_map
        t0 = 0
        bag = Vector{T}()
        for hit in dom_hits
            if hit.t - t0 <= Δt
                push!(bag, hit)
            else
                if length(bag) >= n
                    append!(chits, bag)
                end
                bag = Vector{T}()
            end
            t0 = hit.t
        end
    end
    return chits
end


# Utility
rows(x) = (x[i, :] for i in indices(x,1))


# Math
Base.angle(d1::Direction, d2::Direction) = acos(dot(d1/norm(d1), d2/norm(d2)))
Base.angle(a::T, b::T) where {T<:Union{Hit, PMT, Track}} = Base.angle(a.dir, b.dir)
Base.angle(a::FieldVector{3}, b::Union{Hit, PMT, Track}) = Base.angle(a, b.dir)
Base.angle(a::Union{Hit, PMT, Track}, b::FieldVector{3}) = Base.angle(a.dir, b)

"""
    function pld3(p1, p2, d2)

Calculate the distance between a point (p1) and a line (given by p2 and d2).
"""
function pld3(p1, p2, d2)
    norm(cross(d2, (p2 - p1))) / norm(d2)
end

function matrix(v::Vector)
    m = length(v)
    n = length(v[1])

    M = zeros(m, n)
    i = 1
    for j ∈ 1:n
        for k ∈ 1:m
            M[i] = v[k][j]
            i += 1
        end
    end
    M
end


function svdfit(M)
    com = mean(M, 1)
    subtr = M .- com
    U, S, V = svd(subtr)
    V[:, 1]
end


end # module
