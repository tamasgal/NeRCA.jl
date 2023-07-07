Base.angle(d1, d2) = acos(min(dot(normalize(d1), normalize(d2)), 1))
Base.angle(a::T, b::T) where {T<:Union{KM3io.AbstractCalibratedHit, KM3io.PMT}} = Base.angle(a.dir, b.dir)
Base.angle(a, b::Union{KM3io.AbstractCalibratedHit, KM3io.PMT}) = Base.angle(a, b.dir)
Base.angle(a::Union{KM3io.AbstractCalibratedHit, KM3io.PMT}, b) = Base.angle(a.dir, b)

"""
$(SIGNATURES)

Calculate the cartesian coordinates for given `ϕ`, `θ` and radius `r`.
"""
function cartesian(ϕ, θ; r=1.0)
    return Position(r*sin(θ) * cos(ϕ), r*sin(θ) * sin(ϕ), cos(θ))
end

"""
$(SIGNATURES)

Calculate the distance between a point (p1) and a line (given by p2 and d2).
"""
function pld3(p1, p2, d2)
    norm(cross(d2, (p2 - p1))) / norm(d2)
end


"""
$(SIGNATURES)

Calculate the distance between two lines.
"""
function lld3(P, u, Q, v)
    R = Q - P
    n = cross(u, v)
    return norm(R⋅n) / norm(n)
end


"""
$(SIGNATURES)

Calculate the distance between two tracks.
"""
function lld3(t₁::KM3io.Track, t₂::KM3io.Track)
    return lld3(t₁.pos, t₁.dir, t₂.pos, t₂.dir)   
end


"""
$(SIGNATURES)

Projects a point to a track.
"""
function project(P, t::KM3io.Track)
    A = t.pos
    B = A + t.dir
    project(P, A, B)
end


"""
$(SIGNATURES)

Project P onto a line spanned by A and B.
"""
function project(P, A, B)
    Position(A + ((P-A)⋅(B-A))/((B-A)⋅(B-A)) * (B-A))
end


"""
$(SIGNATURES)

The function to fit time residual plot distributions.
"""
function langauss(x, LA, Lμ, Lσ, GA, Gμ, Gσ, offset)
    LA * pdf(Landau(Lμ, Lσ), x) + GA * pdf(Normal(Gμ, Gσ), x) + offset
end

"""

Creates directions (points on a unit sphere) which are quite evenly distributed
with respect to their space angles using the Fibonacci lattice algorithm. The
axial anisotropy is much smaller compared to a simple latitude-longitude
lattice.

"""
function fibonaccisphere(N)
    ϕ = π * (√5 - 1)  # golden angle in rad

    directions = sizehint!(Direction[], N)
    for i ∈ 0:N-1
        y = 1 - 2(i / (N - 1))
        r = √(1 - y * y)
        θ = ϕ * i
        x = cos(θ) * r
        z = sin(θ) * r
        push!(directions, Direction(x, y, z))
    end
    directions
end
