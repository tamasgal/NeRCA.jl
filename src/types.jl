const ToT = UInt8
const ChannelID = UInt8
const DOMID = UInt32
const Floor = UInt8
const DU = UInt8
const HitTime = Float64

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

struct JMuon
    JENERGY_CHI2::Float64
    JENERGY_ENERGY::Float64
    JGANDALF_BETA0_RAD::Float64
    JGANDALF_BETA1_RAD::Float64
    JGANDALF_CHI2::Float64
    JGANDALF_LAMBDA::Float64
    JGANDALF_NUMBER_OF_HITS::Float64
    JGANDALF_NUMBER_OF_ITERATIONS::Float64
    JMUONENERGY::Bool
    JMUONGANDALF::Bool
    JMUONPREFIT::Bool
    JMUONSIMPLEX::Bool
    JMUONSTART::Bool
    JSTART_LENGTH_METRES::Float64
    JSTART_NPE_MIP::Float64
    JSTART_NPE_MIP_TOTAL::Float64
    dir_x::Float64
    dir_y::Float64
    dir_z::Float64
    energy::Float64
    group_id::Int64
    id::Int64
    length::Float64
    likelihood::Float64
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    rec_type::Int64
    time::Float64
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
    time::Float64
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
    pos::Dict{Int32,Vector{KM3NeT.Position}}
    dir::Dict{Int32,Vector{KM3NeT.Direction}}
    t0::Dict{Int32,Vector{Float64}}
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
    time::HitTime
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
    pos::Position
    dir::Direction
    t0::HitTime
    triggered::Bool
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

Base.isless(lhs::AbstractHit, rhs::AbstractHit) = lhs.time < rhs.time

struct DAQSnapshotHit <: DAQHit
    dom_id::Int32
    channel_id::UInt8
    time::Int32
    tot::UInt8
end

struct DAQTriggeredHit <: DAQHit
    dom_id::Int32
    channel_id::UInt8
    time::Int32
    tot::UInt8
    trigger_mask::Int64
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
    triggered_hits::Vector{DAQTriggeredHit}
    n_snapshot_hits::Int32
    snapshot_hits::Vector{DAQSnapshotHit}
end

Base.show(io::IO, d::DAQEvent) = begin
    print(io, "DAQEvent: $(d.n_triggered_hits) triggered hits, " *
              "$(d.n_snapshot_hits) snapshot hits")
end

function read_io(io::IOBuffer, t::T) where T
    length = read(io, Int32)
    type = read(io, Int32)
    det_id = read(io, Int32)
    run_id = read(io, Int32)
    timeslice_id = read(io, Int32)
    timestamp = read(io, Int32)
    ticks = read(io, Int32)
    trigger_counter = read(io, Int64)
    trigger_mask = read(io, Int64)
    overlays = read(io, Int32)

    n_triggered_hits = read(io, Int32)
    triggered_hits = Vector{DAQTriggeredHit}()
    sizehint!(triggered_hits, n_triggered_hits)
    @inbounds for i ∈ 1:n_triggered_hits
        dom_id = read(io, Int32)
        channel_id = read(io, UInt8)
        time = bswap(read(io, Int32))
        tot = read(io, UInt8)
        trigger_mask = read(io, Int64)
        push!(triggered_hits, DAQTriggeredHit(dom_id, channel_id, time, tot, trigger_mask))
    end

    n_snapshot_hits = read(io, Int32)
    snapshot_hits = Vector{DAQSnapshotHit}()
    sizehint!(snapshot_hits, n_snapshot_hits)
    @inbounds for i ∈ 1:n_snapshot_hits
        dom_id = read(io, Int32)
        channel_id = read(io, UInt8)
        time = bswap(read(io, Int32))
        tot = read(io, UInt8)
        push!(snapshot_hits, DAQSnapshotHit(dom_id, channel_id, time, tot))
    end

    DAQEvent(det_id, run_id, timeslice_id, timestamp, ticks, trigger_counter, trigger_mask, overlays, n_triggered_hits, triggered_hits, n_snapshot_hits, snapshot_hits)
end
