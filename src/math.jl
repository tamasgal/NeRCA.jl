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

Calculate the distance between two lines. `P` is a point on the first line,
`u` is the vector of direction of the same. `Q` and `v` in the same manner for the
second line.

"""
function lld3(P, u, Q, v)
    n = cross(u, v)
    # parallel lines, use point to line distance
    all(iszero, n) && return pld3(P, Q, v)

    R = Q - P
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
function fibonaccisphere(N::Int)
    ϕ = π * (√5 - 1)  # golden angle in rad

    directions = sizehint!(Vector{Direction{Float64}}(), N)
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

Creates directions which with a median angular separation of `α` [deg] using
the Fibonacci lattice.

"""
fibonaccisphere(α::Float64) = fibonaccisphere(Int(ceil((195.39/α)^2)))

"""
Create `S` directions inside a cone with an opening angle of `θ` [deg] which points towards `dir`.
"""
function fibonaccicone(dir::Direction{Float64}, S::Integer, θ::Float64)
    N = Int(ceil(S / sin(deg2rad(θ)/2)^2))
    R = rotation_between(SVector(0.0, 1.0, 0.0), dir)  # the Fibonacci lattice starts spiraling around (0, 1, 0)

    ϕ = π * (√5 - 1)  # golden angle in rad
    directions = sizehint!(Vector{Direction{Float64}}(), S)
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

Creates directions with a median angular separation of `α` [deg] inside a cone with an
opening angle of `θ` [deg] pointing towards `dir`.

"""
function fibonaccicone(dir::Direction{Float64}, α::Float64, θ::Float64)
    α > θ && error("The angular separation needs to be less than the opening angle of the cone.")
    N = (195.39 / α)^2  # total number of points on a full sphere, so that we have α separation
    S = Int(ceil(N * sin(deg2rad(θ)/2)^2))
    fibonaccicone(dir, S, θ)
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

Calculates the maximum angle between the given objects which have a direction.

"""
function spread(objects::Vector{T}) where T
    angles = map(combinations(objects, 2)) do obj_pair
        angle(obj_pair...)
    end
    maximum(angles)
end
