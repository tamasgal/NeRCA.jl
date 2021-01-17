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
    ismxshower, is3dmuon, is3dshower, isnb,
    Calibration, calibrate,
    duhits, domhits, nfoldhits, triggered,
    CHClient, CHTag, CHPrefix, CHMessage, subscribe,
    ztplot,
    savefigs,
    most_frequent,
    SingleDUParams, SingleDURecoParams,
    @ip_str

include("definitions/daqdatatypes.jl")
include("definitions/fitparameters.jl")
include("definitions/reconstruction.jl")
include("definitions/trigger.jl")
include("definitions/w2list_genhen.jl")
include("definitions/w2list_gseagen.jl")

include("constants.jl")
include("types.jl")
include("tools.jl")
include("io/root.jl")
include("calibration.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")

end
