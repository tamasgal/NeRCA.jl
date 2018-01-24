const ToT = UInt8
const ChannelID = UInt8
const DOMID = UInt32
const Floor = UInt8
const DU = UInt8
const HitTime = Float32

struct Position <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end

# Base.*(x::Vector3D, y::Vector3D ) = Vector3D(SVector(y)* SVector(x))

struct Direction <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end

Base.show(io::IO, p::Position) = begin
    s = @sprintf "%.1f %.1f %.1f" p.x p.y p.z
    print(io, s)
end

Base.show(io::IO, d::Direction) = begin
    s = @sprintf "%.1f %.1f %.1f" d.x d.y d.z
    print(io, s)
end

struct EventInfo
    det_id::Int32
    frame_index::UInt32
    livetime_sec::UInt64
    mc_id::Int32
    mc_t::Float64
    n_events_gen::UInt64
    n_files_gen::UInt64
    overlays::UInt32
    trigger_counter::UInt64
    trigger_mask::UInt64
    utc_nanoseconds::UInt64
    utc_seconds::UInt64
    weight_w1::Float64
    weight_w2::Float64
    weight_w3::Float64
    run_id::UInt32
    event_id::UInt32
end

EventInfo(event_info::HDF5.HDF5Compound{17}) = EventInfo(event_info.data...)


# Fit
abstract type AbstractRecoTrack end

struct RecoTrack<:AbstractRecoTrack
    dir::Direction
    pos::Position
    time::Float64
end

struct NoRecoTrack<:AbstractRecoTrack end

# MC
struct Track
    bjorken_y::Float64
    dir::Direction
    pos::Position
    E::Float64
    interaction_channel::UInt32
    is_cc::Bool
    length::Float64
    t::Float64
    particle_type::Int32
end

Track(track::HDF5.HDF5Compound{15}) = begin
    d = track.data
    Track(d[1], Direction(d[2:4]...), Position(d[10:12]...),
          d[5], d[7], d[8], d[9], d[13], d[14])
end

Base.show(io::IO, t::Track) = begin
    E = @sprintf "%0.1f" t.E
    bjorken_y = @sprintf "%0.2f" t.bjorken_y
    print(io, "Track: bjorken_y($(bjorken_y)), t($(t.t)), " *
          "pos($(t.pos)), dir($(t.dir)), E($(E)), type($(t.particle_type))")
end

Base.isless(t1::Track, t2::Track) = t1.E < t2.E


# Hardware
struct PMT
    channel_id::ChannelID
    pos::Position
    dir::Direction
end


struct DOM
    id::UInt32
    floor::Floor
    du::DU
    pmts::Vector{PMT}
end

struct Calibration
    det_id::Int32
    pos::Dict{Int32,Vector{KM3NeT.Position}}
    dir::Dict{Int32,Vector{KM3NeT.Direction}}
    t0::Dict{Int32,Vector{Integer}}
    du::Dict{Int32,DU}
    floor::Dict{Int32,Floor}
end

Base.show(io::IO, c::Calibration) = begin
    print(io, "Calibration data for detector '$(c.det_id)' " *
              "with $(length(c.pos)) modules.")
end

# Signal
abstract type AbstractHit end
abstract type DAQHit<:AbstractHit end

struct Hit <: DAQHit
    channel_id::ChannelID
    dom_id::DOMID
    t::HitTime
    tot::ToT
    triggered::Bool
end

struct McHit <: AbstractHit
    a::Float32
    origin::UInt32
    pmt_id::UInt32
    t::HitTime
end

struct CalibratedHit <: DAQHit
    channel_id::ChannelID
    dom_id::UInt32
    du::DU
    floor::Floor
    t::HitTime
    tot::ToT
    triggered::Bool
    pos::Position
    dir::Direction
    t0::HitTime
end

struct TimesliceHit <: DAQHit
    channel_id::Int8
    dom_id::UInt32
    t::Int32
    tot::Int16
end

Hit(hit::HDF5.HDF5Compound{5}) = begin
    Hit(hit.data...)
end

Base.isless(lhs::AbstractHit, rhs::AbstractHit) = lhs.t < rhs.t



Base.show(io::IO, h::DAQHit) = begin
    print(io, "$(typeof(h)): channel_id($(h.channel_id)), t($(h.t)), " *
          "tot($(h.tot)), dom_id($(h.dom_id)), triggered($(h.triggered))")
end


# Event Loop
struct Event
    id::Integer
    info::EventInfo
    hits::Vector{Hit}
#    mc_hits::Vector{McHit}
    mc_tracks::Vector{Track}
end


Base.show(io::IO, e::Event) = begin
    print(io, "Event $(e.id): $(length(e.hits)) hits, $(length(e.mc_tracks)) MC tracks")
end
