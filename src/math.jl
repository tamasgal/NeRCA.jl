"""
    function angle_between(v1, v2)

Calculates the angle between two vectors
"""
function angle_between(v1, v2)
    _v1 = normalize(v1)
    _v2 = normalize(v2)
    angle = acos(dot(_v1, _v2))
end

function azimuth(d::Direction)
    if d.x >= 0 && d.y >= 0  # Quadrant I
        return atan(d.x/d.y)
    end
    if d.x < 0 && d.y >= 0  # Quadrant II
        return atan(d.x/d.y) + 2π
    end
    return atan(d.x/d.y) + π  # Quadrant III and IV
end
