module NeRCA

using LinearAlgebra
using Statistics
using Printf
using Sockets

using DocStringExtensions

using KM3io
import KM3io: Hit, CalibratedHit, Track
export TriggeredHit

using Parameters
using StaticArrays
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
    RecoTrack, NoRecoTrack,
    OMKey, omkey2domid, floordist,
    duhits, domhits, nfoldhits, triggered,
    CHClient, CHTag, CHPrefix, CHMessage, subscribe,
    ztplot,
    savefigs,
    most_frequent, categorize,
    SingleDUParams, SingleDURecoParams,
    @ip_str,
    initdb, streamds, detx  # db.jl

include("constants.jl")
include("math.jl")
include("hits.jl")
include("mc.jl")
include("fit.jl")
include("controlhost.jl")
include("plot.jl")
include("db.jl")

end
