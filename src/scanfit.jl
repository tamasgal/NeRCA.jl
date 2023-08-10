Base.@kwdef struct MuonScanfitParameters
    tmax::Float64 = 25.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
    tmaxextra::Float64 = 18.0  # [ns]
    nmaxhits::Int = 50  # maximum number of hits to use
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    directions::Vector{Direction{Float64}}
end
MuonScanfit(det::Detector) = MuonScanfit(MuonScanfitParameters(), det, fibonaccisphere(1000))
function Base.show(io::IO, m::MuonScanfit)
   print(io, "$(typeof(m)) in $(length(m.directions)) directions")
end

function (msf::MuonScanfit)(hits::Vector{T}) where T<:KM3io.AbstractHit
    l1builder = L1Builder(L1BuilderParameters(msf.params.tmax, false))
    rhits = l1builder(HitR1, msf.detector, hits)

    sort!(rhits)
    unique!(h->h.dom_id, rhits)

    clique = Clique(Match3B(msf.params.roadwidth, msf.params.tmaxextra))
    clusterize!(rhits, clique)

    candidates = Vector{Tuple{Int, Direction}}()

    # TODO threading is bugged
    # Threads.@threads for dir ∈ msf.directions
    for dir ∈ msf.directions
        rhits_copy = copy(rhits)

        clique1D = Clique(Match1D(msf.params.roadwidth, msf.params.tmaxextra))
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits_copy)
            rhits_copy[idx] = @set rhit.pos = R * rhit.pos
        end

        length(rhits_copy) > msf.params.nmaxhits && resize!(rhits_copy, msf.params.nmaxhits)
        sort!(rhits_copy; by=timetoz)

        clusterize!(rhits_copy, clique1D)

        length(rhits_copy) <= 3 && continue  # TODO 3 comes from the number of parameters, retrieve from Line1Z fitter via type!

        # TODO x-y scane
        push!(candidates, (length(rhits_copy), dir))

    end
    candidates
end

abstract type EstimatorModel end

"""

A straight line parallel to the z-axis.

"""
struct Line1Z <: EstimatorModel
    pos::Position{Float64}
    t::Float64
end
Line1Z() = Line1Z(Position(0.0, 0.0, 0.0), 0.0)
"""

Calculate the Chernkov arrival tive for a given position.

"""
function Base.time(lz::Line1Z, pos::Position)
    v = pos - lz.pos
    R = √(v.x*v.x + v.y*v.y)
    lz.t + (v.z + R * KM3io.Constants.KAPPA_WATER) * KM3io.Constants.C_INVERSE
end


struct Line1ZEstimator
    model::Line1Z
    V::
end
