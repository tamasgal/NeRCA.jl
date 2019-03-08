"""
    function angle_between(v1, v2)

Calculates the angle between two vectors
"""
function angle_between(v1, v2)
    _v1 = normalize(v1)
    _v2 = normalize(v2)
    angle = acos(dot(_v1, _v2))
end
