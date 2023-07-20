module NeRCA

using LinearAlgebra
using Statistics
using Printf

using DocStringExtensions

using KM3io

using Parameters
using StaticArrays
using Rotations
using Setfield
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
    Hit, HitL0, HitL1, HitL2, HitR0, HitR1, HitR2,
    L1Builder, L1BuilderParameters, Match3B,
    RecoTrack, NoRecoTrack, dumandfit,
    MuonScanfit, MuonScanfitParameters,
    duhits, nfoldhits,
    ztplot,
    most_frequent, categorize, modulemap,
    SingleDUParams, SingleDURecoParams,
    initdb, streamds, detx,  # db.jl
    rotator,
    fibonaccisphere,
    # from KM3io
    Detector, Direction, Position



@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """

@template TYPES = """
    $(TYPEDEF)

    $(DOCSTRING)

    # Fields
    $(TYPEDFIELDS)
    """


include("math.jl")
include("hits.jl")
include("mc.jl")
include("scanfit.jl")
include("royfit.jl")
include("plot.jl")
include("db.jl")

end
