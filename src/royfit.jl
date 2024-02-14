mutable struct SingleDUParams
    d::Float64
    t::Float64
    z::Float64
    dz::Float64
    t₀::Float64
end

"""
Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(d_closest, t_closest, z_closest, dir_z, t₀; n=KM3io.Constants.INDEX_OF_REFRACTION_WATER)
    c_ns = KM3io.Constants.c / 1e9
    d_γ(z) = n/√(n^2 - 1) * √(d_closest^2 + (z-z_closest)^2 * (1 - dir_z^2))
    t(z) = (t₀) + 1/c_ns * ((z - z_closest)*dir_z + (n^2 - 1)/n * d_γ(z))
    d_γ, t
end

"""
Returns a function which calculates the arrival time of a Cherenkov photon
at a given position.
"""
function make_cherenkov_calculator(sdp::SingleDUParams; n=KM3io.Constants.INDEX_OF_REFRACTION_WATER)
    c_ns = KM3io.Constants.c / 1e9
    d_γ(z) = n/√(n^2 - 1) * √(sdp.d^2 + (z-sdp.z)^2 * (1 - sdp.dz^2))
    t(z) = (sdp.t₀) + 1/c_ns * ((z - sdp.z)*sdp.dz + (n^2 - 1)/n * d_γ(z))
    d_γ, t
end


