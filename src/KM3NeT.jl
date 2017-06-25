module KM3NeT

using HDF5

export RawHit, TimesliceHit, rows, read_hits


# Basic types
immutable Position{T<:AbstractFloat}
    x::T
    y::T
    z::T
end

immutable Direction{T<:AbstractFloat}
    x::T
    y::T
    z::T
end


# Hardware
immutable PMT
    channel_id::UInt8
    pos::Position
    dir::Direction
end


immutable DOM
    id::UInt32
    floor::UInt8
    line::UInt8
    pmts::Array{PMT}
end


# Signal
abstract type Hit

immutable RawHit <: Hit
    channel_id::UInt8
    dom_id::UInt32
    time::Int32
    tot::UInt8
    triggered::Bool
end

RawHit(hit::HDF5.HDF5Compound{5}) = begin
    RawHit(hit.data...)
end

immutable TimesliceHit <: Hit
    channel_id::Int8
    time::Int32
    tot::Int16
end


Base.show(io::IO, h::RawHit) = begin
    print(io, "$(typeof(h)): channel_id($(h.channel_id)), time($(h.time)), " *
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
    data = h5open(filename, "r") do file
        [RawHit.(read(file, "hits/$i")) for i ∈ event_ids]
    end
    return data
end


# Utility
rows(x) = (x[i, :] for i in indices(x,1))


end # module
