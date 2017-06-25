module KM3NeT

using HDF5

export Track, RawHit, TimesliceHit, rows, read_hits


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


# Hardware
struct PMT
    channel_id::UInt8
    pos::Position
    dir::Direction
end


struct DOM
    id::UInt32
    floor::UInt8
    line::UInt8
    pmts::Array{PMT}
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


# Utility
rows(x) = (x[i, :] for i in indices(x,1))


end # module
