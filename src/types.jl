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
    s = Printf.@sprintf "%.2f %.2f %.2f" p.x p.y p.z
    print(io, s)
end

Base.show(io::IO, d::Direction) = begin
    s = Printf.@sprintf "%.2f %.2f %.2f" d.x d.y d.z
    print(io, s)
end

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

struct TimesliceInfo
    frame_index::UInt32
    slice_id::UInt32
    timestamp::UInt32
    nanoseconds::UInt32
    n_frames::UInt32
    group_id::UInt32
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
