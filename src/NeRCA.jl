module NeRCA

using LinearAlgebra
using Statistics
using Printf
using Sockets

using DocStringExtensions
using Parameters
using StaticArrays
using UnROOT
using PGFPlotsX
using Colors
using RecipesBase
using JuMP
using Ipopt
using IntervalArithmetic
using LandauDistribution
using Distributions: pdf, Normal

using CSV
using DataFrames
using HTTP

import Base: read, +, -, *, getindex, length, eltype

export
    OnlineFile, OfflineFile,
    Position, Direction,
    MCEventInfo, MCTrack, RecoTrack, NoRecoTrack,
    AbstractDAQHit, AbstractMCHit,
    Hit, SnapshotHit, CalibratedHit, MCHit, TimesliceHit,
    DAQEvent,
    ismxshower, is3dmuon, is3dshower, isnb,
    Calibration, calibrate, OMKey, omkey2domid, floordist,
    duhits, domhits, nfoldhits, triggered,
    CHClient, CHTag, CHPrefix, CHMessage, subscribe,
    ztplot,
    savefigs,
    most_frequent, categorize,
    SingleDUParams, SingleDURecoParams,
    @ip_str,
    initdb, streamds, detx  # db.jl

# KM3NeT Dataformat definitions
for inc ∈ readdir(joinpath(@__DIR__, "definitions"), join=true)
    !endswith(inc, ".jl") && continue
    include(inc)
end

include("constants.jl")
include("types.jl")
include("tools.jl")
include("io/detx.jl")
include("io/root.jl")
include("calibration.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")
include("db.jl")

end
