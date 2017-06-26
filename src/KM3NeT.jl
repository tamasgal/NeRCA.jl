module KM3NeT

using HDF5

export Track, RawHit, TimesliceHit, rows, read_hits, Position, Direction


# Basic types
struct Position{T<:AbstractFloat}
    x::T
    y::T
    z::T
end

struct Direction{T<:AbstractFloat}
    x::T
    y::T
    z::T
end

# MC
struct Track
    bjorken_y::Float32
    dir::Direction
    pos::Position
    E::Float64
    interaction_channel::UInt8
    is_cc::Bool
    length::Float64
    t::Int32
    particle_type::Int32
end

Track(track::HDF5.HDF5Compound{15}) = begin
    d = track.data
    Track(d[1], Direction(d[2:4]...), Position(d[10:12]...),
          d[5], d[7], d[8], d[9], d[13], d[14])
end

struct EventInfo
    det_id::Int32
    frame_index::UInt32
    livetime_sec::UInt64
    mc_id::Int32
    mc_t::Float64
    n_events_gen::UInt32
    n_files_gen::UInt32
    overlays::UInt32
    trigger_counter::UInt64
    trigger_mask::UInt64
    utc_nanoseconds::UInt64
    utc_seconds::UInt64
    weight_w1::Float64
    weight_w2::Float64
    weight_w3::Float64
    event_id::UInt32
end

EventInfo(event_info::HDF5.HDF5Compound{16}) = EventInfo(event_info.data...)


# Hardware
struct PMT
    channel_id::UInt8
    pos::Position
    dir::Direction
end


struct DOM
    id::UInt32
    floor::UInt8
    du::UInt8
    pmts::Array{PMT}
end

struct Calibration
    det_id::Int32
    pos::Dict{Int32,Vector{KM3NeT.Position}}
    dir::Dict{Int32,Vector{KM3NeT.Direction}}
    t0::Dict{Int32,Vector{Integer}}
end


# Signal
abstract type Hit end

struct RawHit <: Hit
    channel_id::UInt8
    dom_id::UInt32
    t::Int32
    tot::UInt8
    triggered::Bool
end

struct CalibratedHit <: Hit
    channel_id::UInt8
    dom_id::UInt32
    t::Int32
    tot::UInt8
    triggered::Bool
    pos::Position
    dir::Direction
    t0::Int32
end

RawHit(hit::HDF5.HDF5Compound{5}) = begin
    RawHit(hit.data...)
end

struct TimesliceHit <: Hit
    channel_id::Int8
    time::Int32
    tot::Int16
end


Base.show(io::IO, h::RawHit) = begin
    print(io, "$(typeof(h)): channel_id($(h.channel_id)), time($(h.t)), " *
          "tot($(h.tot)), dom_id($(h.dom_id)), triggered($(h.triggered))")
end


# I/O
function read_hits(filename::AbstractString, event_id::Integer)
    data = h5open(filename, "r") do file
        read(file, "hits/$event_id")
    end
    return RawHit.(data)
end

function read_hits(filename::AbstractString,
                   event_ids::Union{Array{T}, UnitRange{T}}) where {T<:Integer}
    hits = Dict{Int, Array{RawHit, 1}}()
    f = h5open(filename, "r")
    for event_id ∈ event_ids
        hits[event_id] = RawHit.(read(f, "hits/$event_id"))
    end
    close(f)
    return hits
end

function read_tracks(filename::AbstractString)
    tracks = Dict{Int, Array{Track, 1}}()
    f = h5open(filename, "r")
    data = read(f, "mc_tracks")
    for d in data
        event_id = d.data[15]
        if !haskey(tracks, event_id)
            tracks[event_id] = Array{Track, 1}()
        end
        push!(tracks[event_id], Track(d))
    end
    close(f)
    return tracks
end

function read_calibration(filename::AbstractString)
    lines = readlines(filename)

    if 'v' ∈ first(lines)
        det_id, version = map(x->parse(Int,x), split(first(lines), 'v'))
        n_doms = parse(Int, lines[4])
        idx = 5
    else
        det_id, n_doms = map(x->parse(Int,x), split(first(lines)))
        version = 1
        idx = 2
    end

    pos = Dict{Int32,Vector{KM3NeT.Position}}()
    dir = Dict{Int32,Vector{KM3NeT.Direction}}()
    t0s = Dict{Int32,Vector{Int32}}()

    for dom ∈ 1:n_doms
        dom_id, du, floor, n_pmts = map(x->parse(Int,x), split(lines[idx]))
        pos[dom_id] = Vector{KM3NeT.Position}()
        dir[dom_id] = Vector{KM3NeT.Direction}()
        t0s[dom_id] = Vector{Int32}()

        for pmt in 1:n_pmts
            l = split(lines[idx+pmt])
            pmt_id = parse(Int,first(l))
            x, y, z, dx, dy, dz = map(x->parse(Float64, x), l[2:7])
            t0 = parse(Int,first(l[8]))
            push!(pos[dom_id], Position(x, y, z))
            push!(dir[dom_id], Direction(dx, dy, dz))
            push!(t0s[dom_id], t0)
        end
        idx += n_pmts + 1
    end

    Calibration(det_id, pos, dir, t0s)
end

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
        c_hit = CalibratedHit(channel_id, dom_id, t, tot, triggered, pos, dir, t0)
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


function read_event_info(filename::AbstractString)
    event_info = Dict{Int32,EventInfo}()
    entries = h5open(filename) do file
        read(file, "event_info")
    end
    for entry in entries
        e = EventInfo(entry.data...)
        event_info[e.event_id] = e
    end
    event_info
end


# Utility
rows(x) = (x[i, :] for i in indices(x,1))


end # module
