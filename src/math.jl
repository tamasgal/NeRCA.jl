Base.angle(d1, d2) = acos(min(dot(normalize(d1), normalize(d2)), 1))
# TODO: type piracy: this needs to go to KM3io
Base.angle(a::T, b::T) where {T<:Union{KM3io.AbstractCalibratedHit, KM3io.PMT}} = Base.angle(a.dir, b.dir)
Base.angle(a, b::Union{KM3io.AbstractCalibratedHit, KM3io.PMT}) = Base.angle(a, b.dir)
Base.angle(a::Union{KM3io.AbstractCalibratedHit, KM3io.PMT}, b) = Base.angle(a.dir, b)

"""
Calculate the cartesian coordinates for given `ϕ`, `θ` and radius `r`.
"""
function cartesian(ϕ, θ; r=1.0)
    return Position(r*sin(θ) * cos(ϕ), r*sin(θ) * sin(ϕ), cos(θ))
end

"""
Calculate the distance between a point (p1) and a line (given by p2 and d2).
"""
function pld3(p1, p2, d2)
    norm(cross(d2, (p2 - p1))) / norm(d2)
end


"""
Calculate the distance between two lines.
"""
function lld3(P, u, Q, v)
    R = Q - P
    n = cross(u, v)
    return norm(R⋅n) / norm(n)
end


"""
Calculate the distance between two tracks.
"""
function lld3(t₁::KM3io.Track, t₂::KM3io.Track)
    return lld3(t₁.pos, t₁.dir, t₂.pos, t₂.dir)   
end


"""
Projects a point to a track.
"""
function project(P, t::KM3io.Track)
    A = t.pos
    B = A + t.dir
    project(P, A, B)
end


"""
Project P onto a line spanned by A and B.
"""
function project(P, A, B)
    Position(A + ((P-A)⋅(B-A))/((B-A)⋅(B-A)) * (B-A))
end


"""
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

    directions = sizehint!(Vector{Direction}(), N)
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

"""
Create `S` directions inside a cone with an opening angle of `θ` which points towards `dir`.
"""
function fibonaccicone(dir::Direction{Float64}, S::Integer, θ)
    N = Int(ceil(S / sin(θ/2)^2))
    R = rotator(Direction(0.0, 1.0, 0.0), dir)  # the Fibonacci lattice starts spiraling around (0, 1, 0)

    ϕ = π * (√5 - 1)  # golden angle in rad
    directions = sizehint!(Vector{Direction}(), S)
    for i ∈ 0:S-1
        y = 1 - 2(i / (N - 1))
        r = √(1 - y * y)
        θ = ϕ * i
        x = cos(θ) * r
        z = sin(θ) * r
        push!(directions, R * Direction(x, y, z))
    end
    directions
end


"""
Creates directions with a given angular spacing in rad.
"""
omega3d(grid::Float64) = omega3d(Direction(0.0, 0.0, 0.0), (0.0, π), grid)

"""
Creates directions with a principal direction `dir`, a polar angle range `θ` [rad] and
an angular spacing `grid` [rad].
"""
function omega3d(dir::Direction{Float64}, θ::Tuple{Float64, Float64}, grid::Float64)
    directions = Vector{Direction}()

    θ_min = first(θ)
    θ_max = last(θ)

    if (θ_min < 0.0)
      θ_min = 0.0
    end
    if (θ_min > π)
      θ_min = π
    end
    if (θ_max < 0.0)
      θ_max = 0.0
    end
    if (θ_max > π)
      θ_max = π
    end

    error("Not implemented yet")
end


"""
A rotation matrix with a given direction, which will point towards the z-axis when rotated.
"""
function rotator(dir::Direction)
    ct = cos(theta(dir))
    st = sin(theta(dir))
    cp = cos(phi(dir))
    sp = sin(phi(dir))

    RotMatrix(ct*cp, -sp, st*cp, ct*sp, cp, st*sp, -st, 0.0, ct)
end

"""
A rotation matrix which rotates `from_dir` to `to_dir`
"""
function rotator(from_dir::T, to_dir::T) where T<:AbstractVector{Float64}
    v = from_dir × to_dir
    c = from_dir ⋅ to_dir
    u = [   0 -v[3]  v[2];
         v[3]     0 -v[1];
        -v[2]  v[1]     0]
    E = Matrix{Float64}(I, 3, 3)
    (c ≈ -1.0) ? (-E) : (E + u + u * u / (1 + c))
end
