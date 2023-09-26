Base.@kwdef struct MuonScanfitParameters
    tmaxlocal::Float64 = 18.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
    nmaxhits::Int = 50  # maximum number of hits to use
    nfits::Int = 1
    ndirections::Int = 1000  # the number of directions to scan over 4π
    nfinedirections::Int = 500  # number of directions in the cone of the fine scan
    θ::Float64 = 7.0  # opening angle of the fine-scan cone
    σ::Float64 = 5.0  # [ns]
    α₁::Float64 = 1.0  # grid angle of the coarse scan
    α₂::Float64 = 0.5  # grid angle of the fine scan
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
function (msf::MuonScanfit)(hits::Vector{T}) where T<:KM3io.AbstractHit  # 35206 allocations, 10.65 MiB

    rhits = msf.coincidencebuilder(HitR1, msf.detector, hits)  # 1201 allocations, 314.61 KiB

    sort!(rhits)
    unique!(h->h.dom_id, rhits)  # 7 allocations

    clusterize!(rhits, Match3B(msf.params.roadwidth, msf.params.tmaxlocal))  # 3 allocations, 400 bytes (mutates)

    # First round on 4π
    candidates = scanfit(msf.params, rhits, msf.directions)  # 19277 allocations

    isempty(candidates) && return candidates
    sort!(candidates, by=m->m.Q; rev=true)

    # Second round on a directed cone pointing towards the previous best direction
    if msf.params.nfinedirections > 0
        most_likely_dir = first(candidates).dir
        directions = fibonaccicone(most_likely_dir, msf.params.nfinedirections, deg2rad(msf.params.θ))
        candidates = scanfit(msf.params, rhits, directions)  # 5312 allocations, 1000.02 KiB

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
    xytsolvers = Channel{XYTSolver}(Threads.nthreads())
    for _ in Threads.nthreads()
        put!(xytsolvers, XYTSolver(params, 1.0))
    end
    chunk_size = max(1, length(directions) ÷ Threads.nthreads())
    chunks = Iterators.partition(directions, chunk_size)

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            xytsolver = take!(xytsolvers)
            results = map(c -> xytsolver(rhits, c), chunk)
            put!(xytsolvers, xytsolver)
            results
        end
    end
    collect(Iterators.flatten(fetch.(tasks)))
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

    pos = sum(h.pos for h ∈ hits) * W  # 0 allocations
    t = 0.0
    lz = Line1Z(pos, t)  # 0 allocations

    t₀ = sum(time(h) for h ∈ hits) * W * KM3io.Constants.C  # +2000 allocations!

    reset!(est)  # 0 allocations

    y₀ = y₁ = y₂ = 0.0
    hit₀ = first(hits)
    xi = hit₀.pos.x - posx(lz)
    yi = hit₀.pos.y - posy(lz)
    ti = (time(hit₀) * KM3io.Constants.C - t₀ - hit₀.pos.z + posz(lz)) / KM3io.Constants.KAPPA_WATER  # 0 allocations

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


    invert!(est.V, est.MINIMAL_SVD_WEIGHT)  # 9000 allocations

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

    w1 = abs(F.S[1])
    w2 = abs(F.S[2])
    w3 = abs(F.S[3])
    w = max(w1, w2, w3) * precision

    @inbounds for idx in eachindex(F.S)
        F.S[idx] = abs(F.S[idx]) >= w ? 1.0 / F.S[idx] : 0.0
    end

    mul!(V, F.U, diagm(F.S) * F.Vt)
end

struct Variance <: FieldVector{4, Float64}
    x::Float64
    y::Float64
    v::Float64
    w::Float64
end

struct CovMatrix
    α::Float64
    σ::Float64
    M::Matrix{Float64}
    V::Vector{Variance}
    CovMatrix(α, σ, N::Int) = new(α, σ, MMatrix{N, N, Float64, N*N}(undef), MVector{N, Variance}(undef))
end

# TODO: generalise hits parameter
function update!(C::CovMatrix, pos::Position, hits; α=1.0, σ=5.0)
    N = length(hits)

    ta = deg2rad(α)
    ct = cos(ta)
    st = sin(ta)

    for (idx, hit) ∈ enumerate(hits)
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

        C.V[idx] = Variance(x, y, v, w)
    end

    @inbounds for i ∈ 1:N
        @inbounds for j ∈ 1:i
            C.M[i, j] = C.V[i] ⋅ C.V[j]
            C.M[j, i] = C.M[i, j]
        end
        C.M[i, i] = C.V[i] ⋅ C.V[i] + σ^2
    end
    C
end

# TODO: generalise hits parameter
timeresvec(lz::Line1Z, hits) = [time(hit) - time(lz, hit.pos) for hit ∈ hits]


"""
A task worker whichs solves for x, y an t for a given set of hits and a direction.
"""
struct XYTSolver
    hits_buffer::Vector{HitR1}
    covmatrix::CovMatrix
    N::Int
    matcher::Match1D
    # TODO: revise passing params since α is redundant
    function XYTSolver(params::MuonScanfitParameters, α::Float64)
        N = params.nmaxhits
        new(Vector{HitR1}(), CovMatrix(α, params.σ, N), N, Match1D(params.roadwidth, params.tmaxlocal))
    end
end

# TODO: revise if we really need to pass params here
function (s::XYTSolver)(hits::Vector{T}, dir::Direction{Float64}) where T<:AbstractReducedHit
    est = Line1ZEstimator(Line1Z(Position(0, 0, 0), 0))
    χ² = Inf
    R = rotator(dir)
    n_initial_hits = length(hits)
    resize!(s.hits_buffer, n_initial_hits)

    for (idx, hit) ∈ enumerate(hits) # rotate hits
        s.hits_buffer[idx] = @set hit.pos = R * hit.pos
    end

    if n_initial_hits > s.N
        # TODO: review this block, here we may need a partial sort
        sort!(s.hits_buffer; by=timetoz)
        resize!(s.hits_buffer, s.N)
    end

    clusterize!(s.hits_buffer, s.matcher)  # 3053 allocations until here

    NDF = length(s.hits_buffer) - est.NUMBER_OF_PARAMETERS
    N = hitcount(s.hits_buffer)

    length(s.hits_buffer) <= est.NUMBER_OF_PARAMETERS && return MuonScanfitCandidate(Position(0, 0, 0), dir, 0, -Inf, 0)

    sort!(s.hits_buffer)

    try
        estimate!(est, s.hits_buffer)  # +11700 allocations
    catch ex
        # if isa(ex, SingularSVDException)
        # @warn "Singular SVD"
        return MuonScanfitCandidate(Position(0, 0, 0), dir, 0, -Inf, 0)
    end

    # TODO: consider creating a "pos()" getter for everything
    update!(s.covmatrix, est.model.pos, s.hits_buffer)  # 0 allocations
    # TODO: this is really ugly... make update!() return the view itself
    V = view(s.covmatrix.M, 1:length(s.hits_buffer), 1:length(s.hits_buffer))
    Y = timeresvec(est.model, s.hits_buffer)  # about 700 extra allocations here
    V⁻¹ = inv(V)  # +3000 allocations
    χ² = transpose(Y) * V⁻¹ * Y  # +700 allocations
    fit_pos = R \ est.model.pos  # 0 allocations

    MuonScanfitCandidate(fit_pos, dir, est.model.t, quality(χ², N, NDF), NDF)
end
