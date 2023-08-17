Base.@kwdef struct MuonScanfitParameters
    tmaxlocal::Float64 = 18.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
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
        est = Line1ZEstimator(Line1Z(Position(0, 0, 0), 0))

        rhits_copy = copy(rhits)

        clique1D = Clique(Match1D(msf.params.roadwidth, msf.params.tmaxlocal))
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits_copy)
            rhits_copy[idx] = @set rhit.pos = R * rhit.pos
        end

        if length(rhits_copy) > msf.params.nmaxhits
            # TODO: review this block, here we may need a partial sort
            resize!(rhits_copy, msf.params.nmaxhits)
            sort!(rhits_copy; by=timetoz)
        end

        clusterize!(rhits_copy, clique1D)

        length(rhits_copy) < est.NUMBER_OF_PARAMETERS && continue

        try
            estimate!(est, rhits_copy)
        catch ex
            # if isa(ex, SingularSVDException)
            continue
        end




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
distance(l::Line1Z, pos::Position) = √distancesquared(l, pos)
distancesquared(l::Line1Z, pos::Position) = (pos.x - posx(l))^2 + (pos.y - posy(l))^2
posx(l::Line1Z) = l.pos.x
posy(l::Line1Z) = l.pos.y
posz(l::Line1Z) = l.pos.z
posz(l::Line1Z, pos::Position) = l.pos.z - distance(l, pos) / KM3io.Constants.TAN_THETA_C_WATER
"""

Calculate the Chernkov arrival tive for a given position.

"""
function Base.time(lz::Line1Z, pos::Position)
    v = pos - lz.pos
    R = √(v.x*v.x + v.y*v.y)
    lz.t + (v.z + R * KM3io.Constants.KAPPA_WATER) * KM3io.Constants.C_INVERSE
end


struct SingularSVDException <: Exception
    message::String
end

mutable struct Line1ZEstimator
    model::Line1Z
    V::MMatrix{3, 3, Float64, 9}
    NUMBER_OF_PARAMETERS::Int
    MINIMAL_SVD_WEIGHT::Float64
    function Line1ZEstimator(model::Line1Z)
        V = zero(MMatrix{3, 3, Float64, 9})
        new(model, V, 3, 1.0e-4)
    end
end
posx(est::Line1ZEstimator) = posx(est.model)
posy(est::Line1ZEstimator) = posy(est.model)
posz(est::Line1ZEstimator) = posz(est.model)

function reset!(est::Line1ZEstimator)
    est.V .= 0.0
    est
end

# TODO: generalise for "data" using traits
function estimate!(est::Line1ZEstimator, hits)
    N = length(hits)

    N < est.NUMBER_OF_PARAMETERS && error("Not enough data points, $N points, but we require at least $(est.NUMBER_OF_PARAMETERS)")

    W = 1.0 / N

    pos = sum(h.pos for h ∈ hits) * W
    t = 0.0

    t₀ = sum(time(h) for h ∈ hits) * W * KM3io.Constants.C

    reset!(est)

    y₀ = y₁ = y₂ = 0.0
    hit₀ = first(hits)
    xi = hit₀.pos.x - posx(est)
    yi = hit₀.pos.y - posy(est)
    ti = (time(hit₀) * KM3io.Constants.C - t₀ - hit₀.pos.z + posz(est)) / KM3io.Constants.KAPPA_WATER

    # starting from the second hit and including the first in the last iteration
    for idx ∈ 2:N+1
        @inbounds hit = idx > N ? first(hits) : hits[idx]
        xj = hit.pos.x - posx(est)
        yj = hit.pos.y - posy(est)
        tj = (time(hit) * KM3io.Constants.C - t₀ - hit.pos.z + posz(est)) / KM3io.Constants.KAPPA_WATER

        dx = xj - xi
        dy = yj - yi
        dt = ti - tj  # opposite sign

        y = (xj + xi) * dx + (yj + yi) * dy + (tj + ti) * dt

        dx *= 2
        dy *= 2
        dt *= 2

        est.V[1, 2] += dx * dx
        est.V[1, 2] += dx * dy
        est.V[1, 3] += dx * dt
        est.V[2, 2] += dy * dy
        est.V[2, 3] += dy * dt
        est.V[3, 3] += dt * dt

        y₀ += dx * y
        y₁ += dy * y
        y₂ += dt * y

        xi = xj
        yi = yj
        ti = tj
    end

    t₀ = KM3io.Constants.C_INVERSE

    est.V[2, 1] = est.V[1, 2]
    est.V[3, 1] = est.V[1, 3]
    est.V[3, 2] = est.V[2, 3]

    F = svd(est.V)

    abs(F.S[2]) <  est.MINIMAL_SVD_WEIGHT * abs(F.S[1]) && throw(SingularSVDException("$F.S"))

    est.V .= F.U * diagm(1 ./ F.S) * F.Vt

    est.model = Line1Z(
        Position(
            est.V[1, 1] * y₀ + est.V[1, 2] * y₁ + est.V[1, 3] * y₂,
            est.V[2, 1] * y₀ + est.V[2, 2] * y₁ + est.V[2, 3] * y₂,
            pos.z
        ),
        (est.V[3, 1] * y₀ + est.V[3, 2] * y₁ + est.V[3, 3] * y₂) * KM3io.Constants.KAPPA_WATER * KM3io.Constants.C_INVERSE + t₀
    )

    est
end
