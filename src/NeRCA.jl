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
using UnROOT
using PGFPlotsX
using RecipesBase
using Colors
using JuMP
using Ipopt
using IntervalArithmetic

import Base: read, +, -, *

export
    OnlineFile, OfflineFile,
    Position, Direction,
    MCEventInfo, MCTrack, RecoTrack, NoRecoTrack,
    AbstractDAQHit, AbstractMCHit,
    Hit, SnapshotHit, CalibratedHit, MCHit, TimesliceHit,
    DAQEvent,
    is_mxshower, is_3dmuon, is_3dshower,
    Calibration, calibrate,
    duhits, domhits, nfoldhits, triggered,
    CHClient, CHTag, CHPrefix, CHMessage, subscribe,
    ztplot,
    savefigs,
    most_frequent,
    SingleDUParams, SingleDURecoParams,
    @ip_str

include("constants.jl")
include("types.jl")
include("tools.jl")
include("io.jl")
include("calibration.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")

end
