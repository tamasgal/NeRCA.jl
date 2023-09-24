Base.@kwdef struct MuonScanfitParameters
    tmaxlocal::Float64 = 18.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
    nmaxhits::Int = 50  # maximum number of hits to use
    nfits::Int = 1
    ndirections::Int = 1000  # the number of directions to scan over 4π
    nfinedirections::Int = 500  # number of directions in the cone of the fine scan
    θ::Float64 = 7.0  # opening angle of the fine-scan cone
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    directions::Vector{Direction{Float64}}
    coincidencebuilder::L1Builder
    function MuonScanfit(params::MuonScanfitParameters, detector::Detector, directions::Vector{Direction{Float64}})
        coincidencebuilder = L1Builder(L1BuilderParameters(params.tmaxlocal, false))
        new(params, detector, directions, coincidencebuilder)
    end
end
function MuonScanfit(params::MuonScanfitParameters, det::Detector)
    MuonScanfit(params, det, fibonaccisphere(params.ndirections))
end
function MuonScanfit(det::Detector)
    params = MuonScanfitParameters()
    MuonScanfit(params, det, fibonaccisphere(params.ndirections))
end
function Base.show(io::IO, m::MuonScanfit)
   print(io, "$(typeof(m)) in $(length(m.directions)) directions.")
end

"""
Performs a Muon track fit for a given event.
"""
(msf::MuonScanfit)(event::DAQEvent) = msf(event.snapshot_hits)

"""
Performs a Muon track fit for a given set of hits (usually snapshot hits).
"""
function (msf::MuonScanfit)(hits::Vector{T}) where T<:KM3io.AbstractHit
    rhits = msf.coincidencebuilder(HitR1, msf.detector, hits)

    sort!(rhits)
    unique!(h->h.dom_id, rhits)

    clique = Clique(Match3B(msf.params.roadwidth, msf.params.tmaxlocal))
    clusterize!(rhits, clique)

    # First round on 4π
    candidates = scanfit(msf.params, rhits, msf.directions)

    isempty(candidates) && return candidates
    sort!(candidates, by=m->m.Q; rev=true)

    # Second round on a directed cone pointing towards the previous best direction
    if msf.params.nfinedirections > 0
        most_likely_dir = first(candidates).dir
        directions = fibonaccicone(most_likely_dir, msf.params.nfinedirections, deg2rad(msf.params.θ))
        candidates = scanfit(msf.params, rhits, directions)

        isempty(candidates) && return candidates
        sort!(candidates, by=m->m.Q; rev=true)
    end

    candidates[1:msf.params.nfits]
end

"""

Performs the scanfit for each given direction and returns a
`Vector{MuonScanfitCandidate}` with all successful fits. The resulting vector can
be empty if none of the directions had enough hits to perform the algorithm.

"""
function scanfit(params::MuonScanfitParameters, rhits::Vector{T}, directions::Vector{Direction{Float64}}) where T<:AbstractReducedHit
    candidates = Vector{Vector{MuonScanfitCandidate}}()

    for i in 1:Threads.nthreads()
        push!(candidates, MuonScanfitCandidate[])
    end

    Threads.@threads for dir ∈ directions
        est = Line1ZEstimator(Line1Z(Position(0, 0, 0), 0))
        χ² = Inf

        rhits_copy = copy(rhits)

        clique1D = Clique(Match1D(params.roadwidth, params.tmaxlocal))
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits_copy)
            rhits_copy[idx] = @set rhit.pos = R * rhit.pos
        end

        if length(rhits_copy) > params.nmaxhits
            # TODO: review this block, here we may need a partial sort
            resize!(rhits_copy, params.nmaxhits)
            sort!(rhits_copy; by=timetoz)
        end

        # TODO: caveat! mutates both rhits_copy and clique1D, this needs a better interface
        clusterize!(rhits_copy, clique1D)

        NDF = length(rhits_copy) - est.NUMBER_OF_PARAMETERS
        N = hitcount(rhits_copy)

        length(rhits_copy) <= est.NUMBER_OF_PARAMETERS && continue

        sort!(rhits_copy)

        try
            estimate!(est, rhits_copy)
        catch ex
            # if isa(ex, SingularSVDException)
            # @warn "Singular SVD"
            continue
        end

        # TODO: consider creating a "pos()" getter for everything
        # TODO: pass alpha and sigma, like V.set(*this, data.begin(), __end1, gridAngle_deg, sigma_ns);  // JMatrixNZ
        V = covmatrix(est.model.pos, rhits_copy)
        Y = timeresvec(est.model, rhits_copy)
        V⁻¹ = inv(V)
        χ² = transpose(Y) * V⁻¹ * Y
        fit_pos = R \ est.model.pos

        push!(candidates[Threads.threadid()],
              MuonScanfitCandidate(fit_pos, dir, est.model.t, quality(χ², N, NDF), NDF)
        )
    end
    vcat(candidates...)
end

struct MuonScanfitCandidate
    pos::Position{Float64}
    dir::Direction{Float64}
    t::Float64
    Q::Float64
    NDF::Int
end

"""
The quality of the fit, the larger the better, as used in e.g. Jpp.
"""
quality(χ², N, NDF) = N  -  0.25 * χ² / NDF

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
    lz = Line1Z(pos, t)

    t₀ = sum(time(h) for h ∈ hits) * W * KM3io.Constants.C

    reset!(est)

    y₀ = y₁ = y₂ = 0.0
    hit₀ = first(hits)
    xi = hit₀.pos.x - posx(lz)
    yi = hit₀.pos.y - posy(lz)
    ti = (time(hit₀) * KM3io.Constants.C - t₀ - hit₀.pos.z + posz(lz)) / KM3io.Constants.KAPPA_WATER

    # starting from the second hit and including the first in the last iteration
    @inbounds for idx ∈ 2:N+1
        hit = idx > N ? first(hits) : hits[idx]
        xj = hit.pos.x - posx(lz)
        yj = hit.pos.y - posy(lz)
        tj = (time(hit) * KM3io.Constants.C - t₀ - hit.pos.z + posz(lz)) / KM3io.Constants.KAPPA_WATER

        dx = xj - xi
        dy = yj - yi
        dt = ti - tj  # opposite sign

        y = (xj + xi) * dx + (yj + yi) * dy + (tj + ti) * dt

        dx *= 2
        dy *= 2
        dt *= 2

        est.V[1, 1] += dx * dx
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

    t₀ *= KM3io.Constants.C_INVERSE

    @inbounds begin
        est.V[2, 1] = est.V[1, 2]
        est.V[3, 1] = est.V[1, 3]
        est.V[3, 2] = est.V[2, 3]
    end


    invert!(est.V, est.MINIMAL_SVD_WEIGHT)

    @inbounds begin
        est.model = Line1Z(
            Position(
                pos.x + est.V[1, 1] * y₀ + est.V[1, 2] * y₁ + est.V[1, 3] * y₂,
                pos.y + est.V[2, 1] * y₀ + est.V[2, 2] * y₁ + est.V[2, 3] * y₂,
                posz(lz)
            ),
            (est.V[3, 1] * y₀ + est.V[3, 2] * y₁ + est.V[3, 3] * y₂) * KM3io.Constants.KAPPA_WATER * KM3io.Constants.C_INVERSE + t₀
        )
    end

    est
end

"""
Invert matrix in-place with a given precision (clamps eigenvalues to 0 below that).
"""
function invert!(V, precision)
    F = svd(V)

    abs(F.S[2]) <  precision * abs(F.S[1]) && throw(SingularSVDException("$F.S"))

    w = max(map(abs, F.S)...) * precision

    @inbounds for idx in eachindex(F.S)
        F.S[idx] = abs(F.S[idx]) >= w ? 1.0 / F.S[idx] : 0.0
    end

    V .= F.U * diagm(F.S) * F.Vt
    V
end

struct Variance <: FieldVector{4, Float64}
    x::Float64
    y::Float64
    v::Float64
    w::Float64
end

# TODO: generalise hits parameter
function covmatrix(pos::Position, hits; α=1.0, σ=5.0)
    N = length(hits)
    M = Matrix{Float64}(undef, N, N)
    variances = sizehint!(Vector{Variance}(), N)

    ta = deg2rad(α)
    ct = cos(ta)
    st = sin(ta)

    for hit ∈ hits
        dx, dy, dz = hit.pos - pos
        R = √(dx^2 + dy^2)

        x = y = ta * KM3io.Constants.KAPPA_WATER * KM3io.Constants.C_INVERSE
        v = w = ta * KM3io.Constants.C_INVERSE

        if R != 0.0
          x *= dx / R
          y *= dy / R
        end

        x *= (dz * ct - dx * st)
        y *= (dz * ct - dy * st)
        v *= -(dx * ct + dz * st)
        w *= -(dy * ct + dz * st)

        push!(variances, Variance(x, y, v, w))
    end

    @inbounds for i ∈ 1:N
        @inbounds for j ∈ 1:i
            M[i, j] = variances[i] ⋅ variances[j]
            M[j, i] = M[i, j]
        end
        M[i, i] = variances[i] ⋅ variances[i] + σ^2
    end
    M
end

# TODO: generalise hits parameter
timeresvec(lz::Line1Z, hits) = [time(hit) - time(lz, hit.pos) for hit ∈ hits]
