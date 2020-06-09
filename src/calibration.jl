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
