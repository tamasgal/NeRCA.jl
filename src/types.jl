const ToT = UInt8
const ChannelID = UInt8
const DOMID = UInt32
const Floor = UInt8
const DU = UInt8
const HitTime = Float64
const TriggerMask = Int64


struct OMKey
    du::DU
    floor::Floor
end

struct Position{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T

    Position(x::T, y::T, z::T) where {T} = new{T}(x, y, z)
end
Position(x, y, z) = Position(promote(x, y, z)...)

struct Direction{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end
Direction(x, y, z) = Direction(promote(x, y, z)...)
Direction(ϕ, θ) = Direction(cos(ϕ)*cos(θ), sin(ϕ)*cos(θ), sin(θ))

struct MCEventInfo
    event_id::Int64
    group_id::Int64
    mc_id::Int64
    mc_time::Float64
    nanoseconds::Int64
    run_id::Int64
    timestamp::Int64
    weight_w1::Float64
    weight_w2::Float64
    weight_w3::Float64
    weight_w4::Float64
end

Base.show(io::IO, e::MCEventInfo) = begin
    print(io, "MCEventInfo: id($(e.event_id)), mc_id($(e.mc_id)), " *
          "mc_time($(e.mc_time) ($(e.nanoseconds)), ts($(e.timestamp)), run_id($(e.run_id))")
end

struct DAQEventInfo
    det_id::Int64
    event_id::Int64
    frame_index::Int64
    group_id::Int64
    mc_run_id::Int64
    mc_time::Float64
    nanoseconds::Int64
    overlays::Int64
    run_id::Int64
    timestamp::Int64
    trigger_counter::Int64
    trigger_mask::Int64
    weight_w1::Float64
    weight_w2::Float64
    weight_w3::Float64
    weight_w4::Float64
end

Base.show(io::IO, e::DAQEventInfo) = begin
    print(io, "DAQEventInfo: id($(e.event_id)), det_id($(e.det_id)), " *
          "timestamp($(e.timestamp) ($(e.nanoseconds)), run_id($(e.run_id))")
end


struct TimesliceInfo
    frame_index::UInt32
    slice_id::UInt32
    timestamp::UInt32
    nanoseconds::UInt32
    n_frames::UInt32
    group_id::UInt32
end

struct Track
    dir::Direction
    pos::Position
    time
end

# Fit
abstract type AbstractRecoTrack end

struct RecoTrack<:AbstractRecoTrack
    dir::Direction
    pos::Position
    time::Float64
end

struct NoRecoTrack<:AbstractRecoTrack end

# MC
struct MCTrack
    bjorken_y::Float64
    dir_x::Float64
    dir_y::Float64
    dir_z::Float64
    E::Float64
    group_id::Int64
    id::Int64
    interaction_channel::Int64
    is_cc::Bool
    length::Float64
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    t::Float64
    particle_type::Int64
end

Track(t::MCTrack) = Track([t.dir_x, t.dir_y, t.dir_z], [t.pos_x, t.pos_y, t.pos_z], t.t)


Base.show(io::IO, t::MCTrack) = begin
    E = Printf.@sprintf "%0.1f" t.E
    bjorken_y = Printf.@sprintf "%0.2f" t.bjorken_y
    print(io, "MCTrack: bjorken_y($(bjorken_y)), t($(t.t)), " *
              "pos($(t.pos_x), $(t.pos_y), $(t.pos_z)), " *
              "dir($(t.dir_x), $(t.dir_y), $(t.dir_z)), " *
              "E($(E)), type($(t.particle_type))")
end


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
    pos::Dict{Int32,Vector{NeRCA.Position}}
    dir::Dict{Int32,Vector{NeRCA.Direction}}
    t0::Dict{Int32,Vector{Float64}}
    du::Dict{Int32,DU}
    floor::Dict{Int32,Floor}
    omkeys::Dict{OMKey,Int32}  # maps OMKey to Module ID
    pmts::Dict{Int32,Tuple{Int32,UInt8}}
    max_z
    n_dus
end

Base.show(io::IO, c::Calibration) = begin
    print(io, "Calibration data for detector '$(c.det_id)' " *
              "with $(length(c.pos)) modules.")
end

# Signal
abstract type AbstractHit end
abstract type AbstractDAQHit<:AbstractHit end
abstract type AbstractMCHit<:AbstractHit end


Base.isless(lhs::AbstractHit, rhs::AbstractHit) = lhs.t < rhs.t


struct Hit <: AbstractDAQHit
    channel_id::ChannelID
    dom_id::DOMID
    t::HitTime
    tot::ToT
    trigger_mask::TriggerMask
end


struct SnapshotHit <: AbstractDAQHit
    channel_id::ChannelID
    dom_id::DOMID
    t::HitTime
    tot::ToT
end

struct MCHit <: AbstractHit
    a::Float32
    origin::UInt32
    pmt_id::UInt32
    t::HitTime
end

mutable struct Multiplicity
    count::Int32
    id::Int64
end

struct CalibratedHit <: AbstractDAQHit
    channel_id::ChannelID
    dom_id::UInt32
    du::DU
    floor::Floor
    t::HitTime
    tot::ToT
    pos::Position
    dir::Direction
    t0::HitTime
    trigger_mask::TriggerMask
    multiplicity::Multiplicity
end


struct TimesliceHit <: AbstractDAQHit
    channel_id::Int8
    dom_id::UInt32
    t::Int32
    tot::Int16
end


struct TriggeredHit <: AbstractDAQHit
    dom_id::Int32
    channel_id::UInt8
    t::Int32
    tot::UInt8
    trigger_mask::TriggerMask
end

struct DAQEvent
    det_id::Int32
    run_id::Int32
    timeslice_id::Int32
    timestamp::Int32
    ticks::Int32
    trigger_counter::Int64
    trigger_mask::Int64
    overlays::Int32
    n_triggered_hits::Int32
    triggered_hits::Vector{TriggeredHit}
    n_hits::Int32
    hits::Vector{Hit}
end

Base.show(io::IO, d::DAQEvent) = begin
    print(io, "DAQEvent: $(d.n_triggered_hits) triggered hits, " *
              "$(d.n_hits) snapshot hits")
end

