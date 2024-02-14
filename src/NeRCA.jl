module NeRCA

using LinearAlgebra
using Statistics
using Printf

using DocStringExtensions

using KM3io

using StaticArrays
using Rotations
using Combinatorics
using Setfield
using Colors
using LandauDistribution
using Distributions: pdf, Normal

using CSV
using DataFrames
using HTTP

import Base: read, +, -, *, getindex, length, eltype

export
    Hit, HitL0, HitL1, HitL2, HitR0, HitR1, HitR2,
    L1Builder, L1BuilderParameters, Match3B, Match1D,
    Line1Z, Line1ZEstimator,
    dumandfit,
    MuonScanfit, MuonScanfitCandidate, MuonScanfitParameters, timetoz,
    duhits, nfoldhits,
    most_frequent, categorize, modulemap,
    initdb, streamds, detx,  # db.jl
    fibonaccisphere, fibonaccicone, rotator, spread,
    # re-export from KM3io
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
include("dumandfit.jl")
include("royfit.jl")
include("db.jl")

end
