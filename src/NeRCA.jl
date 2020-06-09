module NeRCA

@time using LinearAlgebra
@time using Statistics
@time using Printf
@time using Sockets

@time using DocStringExtensions
@time using DrWatson
@time using Parameters
@time using DataStructures
@time using StaticArrays
@time using HDF5
@time using UnROOT
@time using PGFPlotsX
@time using RecipesBase
@time using PlotThemes
@time using JuMP
@time using Ipopt
@time using IntervalArithmetic

@time import Base: +, -, *

println("Done loading libraries...")

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
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")

end
