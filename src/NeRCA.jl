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
    OnlineFile, OfflineFile,
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
include("calibration.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")

end
