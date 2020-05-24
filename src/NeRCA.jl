module NeRCA

using LinearAlgebra
using Statistics
using Printf
using Sockets

using DocStringExtensions
using DrWatson
using Parameters
using DataStructures
using StaticArrays
using HDF5
using UnROOT
using PGFPlotsX
using RecipesBase
using PlotThemes
using JuMP
using Ipopt
using IntervalArithmetic

import Base: +, -, *

export
    Position, Direction,
    MCEventInfo, MCTrack, RecoTrack, NoRecoTrack,
    AbstractDAQHit, AbstractMCHit,
    Hit, SnapshotHit, CalibratedHit, MCHit, TimesliceHit,
    DAQEvent, read_io,
    is_mxshower, is_3dmuon, is_3dshower,
    read_indices, read_hits, read_tracks, read_mctracks,
    read_calibration, read_event_info,
    Calibration, calibrate,
    duhits, domhits, nfoldhits, triggered,
    CHClient, CHTag, CHPrefix, CHMessage, subscribe,
    ztplot,
    savefigs,
    most_frequent,
    @ip_str


include("constants.jl")
include("types.jl")
include("tools.jl")
include("io.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")


"""
$(SIGNATURES)

Apply geometry and time calibration to given hits.
"""
function calibrate(calibration::Calibration, hits)
    calibrated_hits = Vector{CalibratedHit}()
    has_triggermask = hasfield(eltype(hits), :trigger_mask)
    for hit in hits
        dom_id = hit.dom_id
        channel_id = hit.channel_id
        tot = hit.tot
        pos = calibration.pos[dom_id][channel_id+1]
        dir = calibration.dir[dom_id][channel_id+1]
        t0 = calibration.t0[dom_id][channel_id+1]
        t = hit.t + t0
        du = calibration.du[dom_id]
        floor = calibration.floor[dom_id]
        trigger_mask = 0
        if has_triggermask
            trigger_mask = hit.trigger_mask
        end
        c_hit = CalibratedHit(channel_id, dom_id, du, floor, t, tot,
                              pos, dir, t0, trigger_mask, Multiplicity(0,0))
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


"""
$(SIGNATURES)

Apply geometry and time calibration to given mc_hits.
"""
function calibrate(mc_hits::Vector{T},
                   calibration::Calibration,
                   event_info::Union{Nothing,MCEventInfo,DAQEventInfo} = nothing) where {T<:MCHit}
    calibrated_hits = Vector{CalibratedHit}()
    if event_info != nothing
        mctime = make_mc_time_converter(event_info)
    else
        mctime = x->x
    end
    for hit in mc_hits
        omkey = calibration.omkeys[hit.pmt_id]
        dom_id = omkey.dom_id
        channel_id = omkey.channel_id
        tot = hit.a
        pos = calibration.pos[dom_id][channel_id+1]
        dir = calibration.dir[dom_id][channel_id+1]
        t0 = calibration.t0[dom_id][channel_id+1]
        t = mctime(hit.t)
        du = calibration.du[dom_id]
        floor = calibration.floor[dom_id]
        triggered = false
        if T === SnapshotHit
            triggered = hit.triggered
        end
        if T === TriggeredHit
            triggered = hit.trigger_mask > 0
        end
        c_hit = CalibratedHit(channel_id, dom_id, du, floor, t, tot,
                              pos, dir, t0, triggered, Multiplicity(0,0))
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


# Math
Base.angle(d1::Direction, d2::Direction) = acos(dot(d1/norm(d1), d2/norm(d2)))
Base.angle(a::T, b::T) where {T<:Union{CalibratedHit, PMT, MCTrack}} = Base.angle(a.dir, b.dir)
Base.angle(a::FieldVector{3}, b::Union{CalibratedHit, PMT, MCTrack}) = Base.angle(a, b.dir)
Base.angle(a::Union{CalibratedHit, PMT, MCTrack}, b::FieldVector{3}) = Base.angle(a.dir, b)


end
