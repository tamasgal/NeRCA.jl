"""
    function angle_between(v1, v2)

Calculates the angle between two vectors
"""
function angle_between(v1, v2)
    _v1 = normalize(v1)
    _v2 = normalize(v2)
    angle = acos(dot(_v1, _v2))
end


"""
    function azimuth(d::Direction)

Calculate the azimuth angle for a given direction.
"""
function azimuth(d::Direction)
    if d.x >= 0 && d.y >= 0  # Quadrant I
        return atan(d.x/d.y)
    end
    if d.x < 0 && d.y >= 0  # Quadrant II
        return atan(d.x/d.y) + 2π
    end
    return atan(d.x/d.y) + π  # Quadrant III and IV
end


"""
    function lld3(P, u, Q, v)

Calculate the distance between two lines.
"""
function lld3(P, u, Q, v)
    R = Q - P
    n = cross(u, v)
    return norm(R⋅n) / norm(n)
end


"""
    function lld3(t₁::T, t₂::T) where T::Track

Calculate the distance between two tracks.
"""
function lld3(t₁::Track, t₂::Track)
    return lld3(t₁.pos, t₁.dir, t₂.pos, t₂.dir)   
end

