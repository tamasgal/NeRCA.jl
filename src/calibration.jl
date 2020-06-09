"""
$(SIGNATURES)

Apply geometry and time calibration to given hits.
"""
function calibrate(calibration::Calibration, hits)
    calibrated_hits = Vector{CalibratedHit}()
    has_trigger_mask = :trigger_mask ∈ fieldnames(eltype(hits))
    for hit in hits
        dom_id = hit.dom_id
        channel_id = hit.channel_id
        tot = hit.tot
        pos = calibration.pos[dom_id][channel_id+1]
        dir = calibration.dir[dom_id][channel_id+1]
        t0 = calibration.t0[dom_id][channel_id+1]
        t = hit.t + t0
        du = calibration.du[dom_id]
        floor = calibration.floor[dom_id]
        trigger_mask = has_trigger_mask ? hit.trigger_mask : 0
        c_hit = CalibratedHit(channel_id, dom_id, du, floor, t, tot,
                              pos, dir, t0, trigger_mask, Multiplicity(0,0))
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


"""
$(SIGNATURES)

Apply geometry and time calibration to given mc_hits.
"""
function calibrate(mc_hits::Vector{T},
                   calibration::Calibration,
                   event_info::Union{Nothing,MCEventInfo,DAQEventInfo} = nothing) where {T<:MCHit}
    calibrated_hits = Vector{CalibratedHit}()
    if event_info != nothing
        mctime = make_mc_time_converter(event_info)
    else
        mctime = x->x
    end
    for hit in mc_hits
        omkey = calibration.omkeys[hit.pmt_id]
        dom_id = omkey.dom_id
        channel_id = omkey.channel_id
        tot = hit.a
        pos = calibration.pos[dom_id][channel_id+1]
        dir = calibration.dir[dom_id][channel_id+1]
        t0 = calibration.t0[dom_id][channel_id+1]
        t = mctime(hit.t)
        du = calibration.du[dom_id]
        floor = calibration.floor[dom_id]
        c_hit = CalibratedHit(channel_id, dom_id, du, floor, t, tot,
                              pos, dir, t0, hit.triggered, Multiplicity(0,0))
        push!(calibrated_hits, c_hit)
    end
    calibrated_hits
end


const SLEWS = [
    8.01, 7.52, 7.05, 6.59, 6.15, 5.74, 5.33, 4.95, 4.58, 4.22, 3.89, 3.56,
    3.25, 2.95, 2.66, 2.39, 2.12, 1.87, 1.63, 1.40, 1.19, 0.98, 0.78, 0.60,
    0.41, 0.24, 0.07, -0.10, -0.27, -0.43, -0.59, -0.75, -0.91, -1.08,
    -1.24, -1.41, -1.56, -1.71, -1.85, -1.98, -2.11, -2.23, -2.35, -2.47,
    -2.58, -2.69, -2.79, -2.89, -2.99, -3.09, -3.19, -3.28, -3.37, -3.46,
    -3.55, -3.64, -3.72, -3.80, -3.88, -3.96, -4.04, -4.12, -4.20, -4.27,
    -4.35, -4.42, -4.49, -4.56, -4.63, -4.70, -4.77, -4.84, -4.90, -4.97,
    -5.03, -5.10, -5.16, -5.22, -5.28, -5.34, -5.40, -5.46, -5.52, -5.58,
    -5.63, -5.69, -5.74, -5.80, -5.85, -5.91, -5.96, -6.01, -6.06, -6.11,
    -6.16, -6.21, -6.26, -6.31, -6.36, -6.41, -6.45, -6.50, -6.55, -6.59,
    -6.64, -6.68, -6.72, -6.77, -6.81, -6.85, -6.89, -6.93, -6.98, -7.02,
    -7.06, -7.09, -7.13, -7.17, -7.21, -7.25, -7.28, -7.32, -7.36, -7.39,
    -7.43, -7.46, -7.50, -7.53, -7.57, -7.60, -7.63, -7.66, -7.70, -7.73,
    -7.76, -7.79, -7.82, -7.85, -7.88, -7.91, -7.94, -7.97, -7.99, -8.02,
    -8.05, -8.07, -8.10, -8.13, -8.15, -8.18, -8.20, -8.23, -8.25, -8.28,
    -8.30, -8.32, -8.34, -8.37, -8.39, -8.41, -8.43, -8.45, -8.47, -8.49,
    -8.51, -8.53, -8.55, -8.57, -8.59, -8.61, -8.62, -8.64, -8.66, -8.67,
    -8.69, -8.70, -8.72, -8.74, -8.75, -8.76, -8.78, -8.79, -8.81, -8.82,
    -8.83, -8.84, -8.86, -8.87, -8.88, -8.89, -8.90, -8.92, -8.93, -8.94,
    -8.95, -8.96, -8.97, -8.98, -9.00, -9.01, -9.02, -9.04, -9.04, -9.04,
    -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04,
    -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04,
    -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04,
    -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04,
    -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04, -9.04,
    -9.04, -9.04
]


"""
$(SIGNATURES)

Return the time slewing for a ToT.
"""
slew(tot) = SLEWS[tot + 1]
slew(hit::AbstractHit) = slew(hit.tot)
