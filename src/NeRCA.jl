module NeRCA

using LinearAlgebra
using Statistics
using Printf

using DocStringExtensions

using KM3io

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
    RecoTrack, NoRecoTrack, dumandfit,
    duhits, domhits, nfoldhits,
    ztplot,
    most_frequent, categorize,
    SingleDUParams, SingleDURecoParams,
    initdb, streamds, detx  # db.jl

include("math.jl")
include("hits.jl")
include("mc.jl")
include("prefit.jl")
include("fit.jl")
include("plot.jl")
include("db.jl")

end
