Base.@kwdef struct MuonScanfitParameters
    tmaxlocal::Float64 = 18.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
    nmaxhits::Int = 50  # maximum number of hits to use
    nfits::Int = 1
    nprefits::Int = 10
    σ::Float64 = 5.0  # [ns]
    α₁::Float64 = 7.0  # grid angle of the coarse scan
    α₂::Float64 = 0.5  # grid angle of the fine scan
    θ::Float64 = 3.5  # opening angle of the fine-scan cone
end


"""

A container of directions with additionial information about their median
angular separation.

"""
struct DirectionSet
    directions::Vector{Direction{Float64}}
    angular_separation::Float64
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    coarsedirections::DirectionSet
    coincidencebuilder::L1Builder
    function MuonScanfit(params::MuonScanfitParameters, detector::Detector)
        coincidencebuilder = L1Builder(L1BuilderParameters(params.tmaxlocal, false))
        new(params, detector, DirectionSet(fibonaccisphere(params.α₁), params.α₁), coincidencebuilder)
    end
end
MuonScanfit(det::Detector) = MuonScanfit(MuonScanfitParameters(), det)
function Base.show(io::IO, m::MuonScanfit)
    print(io, "$(typeof(m)) with a coarse scan of $(m.params.α₁)ᵒ and a fine scan of $(m.params.α₂)ᵒ.")
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

    clusterize!(rhits, Match3B(msf.params.roadwidth, msf.params.tmaxlocal))

    # First round on 4π
    candidates = scanfit(msf.params, rhits, msf.coarsedirections)
    isempty(candidates) && return candidates
    sort!(candidates, by=m->m.Q; rev=true)

    # Second round on directed cones pointing towards the previous best directions
    # TODO: currently disabled until all the allocations are minimised
    # here, reusing a Vector{Direction} (attached to msf as buffer) might be a good idea.
    # By doing so, we need `fibonaccicone!` and `fibonaccisphere!` as mutating functions
    if msf.params.nprefits > 0
        directions = Vector{Vector{Direction{Float64}}}()
        for idx in 1:min(msf.params.nprefits, length(candidates))
            most_likely_dir = candidates[idx].dir
            push!(directions, fibonaccicone(most_likely_dir, msf.params.α₂, msf.params.θ))
        end
        directionset = DirectionSet(vcat(directions...), msf.params.α₂)
        candidates = scanfit(msf.params, rhits, directionset)

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
function scanfit(params::MuonScanfitParameters, rhits::Vector{T}, directionset::DirectionSet) where T<:AbstractReducedHit
    chunk_size = max(1, length(directionset.directions) ÷ Threads.nthreads())
    chunks = Iterators.partition(directionset.directions, chunk_size)

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            xytsolver = XYTSolver(params.nmaxhits, params.roadwidth, params.tmaxlocal, params.σ)
            results = [xytsolver(rhits, c, directionset.angular_separation) for c in chunk]
            results
        end
    end

    mapreduce(fetch, vcat, tasks)
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

    yvec = zeros(MVector{3})
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

        yvec[1] += dx * y
        yvec[2] += dy * y
        yvec[3] += dt * y

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


    # Hermitian is needed for typestability!
    wvec, evecs = invert2!(Hermitian(est.V), est.MINIMAL_SVD_WEIGHT)
    yvec2 = (evecs' * yvec)
    yvec2 .*= wvec
    mul!(yvec, evecs, yvec2)
    #yvec = evecs * (diagm(wvec) * (evecs' * yvec))

    @inbounds begin
        est.model = Line1Z(
            Position(
                pos.x + yvec[1],
                pos.y + yvec[2],
                posz(lz)
            ),
            yvec[3] * KM3io.Constants.KAPPA_WATER * KM3io.Constants.C_INVERSE + t₀
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

@inline function invert2!(V, precision)
    evals, evecs = eigen(V)

    abs(evals[2]) <  precision * abs(evals[1]) && throw(SingularSVDException("$evals"))

    w = maximum(abs, evals) * precision

    wvec = ifelse.(abs.(evals) .>= w, inv.(evals), zero(float(eltype(evals))))

    return wvec, evecs
    #mul!(V, evecs, diagm(wvec) * evecs')
end

struct Variance <: FieldVector{4, Float64}
    x::Float64
    y::Float64
    v::Float64
    w::Float64
end

struct CovMatrix
    M::Matrix{Float64}
    V::Vector{Variance}
    σ::Float64
    CovMatrix(N::Int, σ::Float64) = new(MMatrix{N, N, Float64, N*N}(undef), MVector{N, Variance}(undef), σ)
end

# TODO: generalise hits parameter
function update!(C::CovMatrix, pos::Position, hits, α::Float64)
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
        C.M[i, i] = C.V[i] ⋅ C.V[i] + C.σ^2
    end
    C
end

# TODO: generalise hits parameter
timeresvec(lz::Line1Z, hits) = [time(hit) - time(lz, hit.pos) for hit ∈ hits]
function timeresvec!(v::AbstractArray{Float64}, lz::Line1Z, hits)
    for (idx, hit) ∈ enumerate(hits)
        v[idx] = time(hit) - time(lz, hit.pos)
    end
    v
end


"""
A task worker whichs solves for x, y an t for a given set of hits and a direction.
"""
struct XYTSolver
    hits_buffer::Vector{HitR1}
    covmatrix::CovMatrix
    timeresvec::Vector{Float64}
    nmaxhits::Int
    matcher::Match1D
    est::Line1ZEstimator

    function XYTSolver(nmaxhits::Int, roadwidth::Float64, tmaxlocal::Float64, σ::Float64)
        new(Vector{HitR1}(), CovMatrix(nmaxhits, σ), Vector{Float64}(), nmaxhits, Match1D(roadwidth, tmaxlocal),
            Line1ZEstimator(Line1Z(Position(0, 0, 0), 0))
        )
    end
end

function (s::XYTSolver)(hits::Vector{T}, dir::Direction{Float64}, α::Float64) where T<:AbstractReducedHit
    χ² = Inf
    R = rotator(dir)
    n_initial_hits = length(hits)
    resize!(s.hits_buffer, n_initial_hits)

    for (idx, hit) ∈ enumerate(hits) # rotate hits
        s.hits_buffer[idx] = @set hit.pos = R * hit.pos
    end

    if n_initial_hits > s.nmaxhits
        sort!(s.hits_buffer; by=timetoz, alg=PartialQuickSort(s.nmaxhits))
        resize!(s.hits_buffer, s.nmaxhits)
    end

    clusterize!(s.hits_buffer, s.matcher)

    hits = s.hits_buffer  # just for convenience
    n_final_hits = length(hits)

    n_final_hits <= s.est.NUMBER_OF_PARAMETERS && return MuonScanfitCandidate(Position(0, 0, 0), dir, 0, -Inf, 0)

    NDF = n_final_hits - s.est.NUMBER_OF_PARAMETERS
    N = hitcount(hits)
    sort!(hits)

    try
        estimate!(s.est, hits)
    catch ex
        # isa(ex, SingularSVDException) && @warn "Singular SVD"
        return MuonScanfitCandidate(Position(0, 0, 0), dir, 0, -Inf, 0)
    end

    # TODO: consider creating a "pos()" getter for everything
    update!(s.covmatrix, s.est.model.pos, hits, α)
    # TODO: this is really ugly... make update!() return the view itself maybe?
    V = view(s.covmatrix.M, 1:n_final_hits, 1:n_final_hits)

    n_final_hits > length(s.timeresvec) && resize!(s.timeresvec, n_final_hits)
    # TODO: better name for this function
    timeresvec!(s.timeresvec, s.est.model, hits)

    Y = view(s.timeresvec, 1:n_final_hits)  # only take the relevant part of the buffer
    χ² = dot(Y, V \ Y)
    fit_pos = R \ s.est.model.pos

    MuonScanfitCandidate(fit_pos, dir, s.est.model.t, quality(χ², N, NDF), NDF)
end
