"""
$(SIGNATURES)

Combine snapshot and triggered hits to a single hits-vector.

This should be used to transfer the trigger information to the
snapshot hits from a DAQEvent. The triggered hits are a subset
of the snapshot hits.

"""
function combine(snapshot_hits::Vector{KM3io.SnapshotHit}, triggered_hits::Vector{KM3io.TriggeredHit})
    triggermasks = Dict{Tuple{UInt8, Int32, Int32, UInt8}, Int64}()
    for hit ∈ triggered_hits
        triggermasks[(hit.channel_id, hit.dom_id, hit.t, hit.tot)] = hit.trigger_mask
    end
    n = length(snapshot_hits)
    hits = sizehint!(Vector{TriggeredHit}(), n)
    for hit in snapshot_hits
        channel_id = hit.channel_id
        dom_id = hit.dom_id
        t = hit.t
        tot = hit.tot
        triggermask = get(triggermasks, (channel_id, dom_id, t, tot), 0)
        push!(hits, TriggeredHit(dom_id, channel_id, t, tot, triggermask))
    end
    hits
end


"""
$(SIGNATURES)

Create a `Vector` with hits contributing to `n`-fold coincidences within a time
window of Δt.
"""
function nfoldhits(hits::Vector{T}, Δt, n) where {T<:KM3io.AbstractDAQHit}
    hit_map = domhits(hits)
    chits = Vector{T}()
    for (dom_id, dom_hits) ∈ hit_map
        bag = Vector{T}()
        push!(bag, dom_hits[1])
        t0 = dom_hits[1].t
        for hit in dom_hits[2:end]
            if hit.t - t0 > Δt
                if length(bag) >= n
                    append!(chits, bag)
                end
                bag = Vector{T}()
            end
            push!(bag, hit)
            t0 = hit.t
        end
    end
    return chits
end


"""
$(SIGNATURES)

Calculate the multiplicities for a given time window. Two arrays are
are returned, one contains the multiplicities, the second one the IDs
of the coincidence groups.
The hits should be sorted by time and then by dom_id.
"""
function count_multiplicities(hits::Vector{T}, tmax=20) where {T<:KM3io.AbstractHit}
    n = length(hits)
    mtp = ones(Int32, n)
    cid = zeros(Int32, n)
    idx0 = 1
    _mtp = 1
    _cid = idx0
    t0 = hits[idx0].t
    dom_id = hits[idx0].dom_id
    for i in 2:n
        hit = hits[i]
        if hit.dom_id != dom_id
            mtp[idx0:i-1] .= _mtp
            cid[idx0:i-1] .= _cid
            dom_id = hit.dom_id
            t0 = hit.t
            _mtp = 1
            _cid += 1
            idx0 = i
            continue
        end
        Δt = hit.t - t0
        if Δt > tmax
            mtp[idx0:i] .= _mtp
            cid[idx0:i] .= _cid
            _mtp = 0
            _cid += 1
            idx0 = i
            t0 = hit.t
        end
        _mtp += 1
    end
    mtp[idx0:end] .= _mtp
    cid[idx0:end] .= _cid
    mtp, cid
end

"""
$(SIGNATURES)

Counts the multiplicities and modifies the .multiplicity field of the hits.
Important: the hits have to be sorted by time and then by DOM ID first.
"""
function count_multiplicities!(hits::Vector{KM3io.XCalibratedHit}, tmax=20)
    _mtp = 0
    _cid = 0
    t0 = 0
    dom_id = 0
    hit_buffer = Vector{XCalibratedHit}()

    function process_buffer()
        while !isempty(hit_buffer)
            _hit = pop!(hit_buffer)
            _hit.multiplicity.count = _mtp
            _hit.multiplicity.id = _cid
        end
    end

    function reset()
        _mtp = 1
        _cid += 1
    end

    for hit ∈ hits
        if length(hit_buffer) == 0
            reset()
            push!(hit_buffer, hit)
            t0 = hit.t
            dom_id = hit.dom_id
            continue
        end
        if hit.dom_id != dom_id
            process_buffer()
            push!(hit_buffer, hit)
            t0 = hit.t
            dom_id = hit.dom_id
            reset()
            continue
        end
        Δt = hit.t - t0
        if Δt > tmax
            process_buffer()
            push!(hit_buffer, hit)
            t0 = hit.t
            reset()
        else
            push!(hit_buffer, hit)
            _mtp += 1
        end
    end
    if length(hit_buffer) > 0
        process_buffer()
    end
    return
end


"""
$(SIGNATURES)

Categorise hits by DOM ID and put them into a dictionary of DOM ID=>Vector{Hit}.
"""
@inline domhits(hits) = categorize(:dom_id, hits)


"""
$(SIGNATURES)

Categorise hits by DU and put them into a dictionary of DU=>Vector{Hit}.
"""
@inline duhits(hits::Vector{T}) where {T<:KM3io.XCalibratedHit} = categorize(:du, hits)


"""
$(SIGNATURES)

Return a vector of hits with ToT >= `tot`.
"""
function totcut(hits::Vector{T}, tot) where {T<:KM3io.AbstractDAQHit}
    return filter(h->h.tot >= tot, hits)
end


"""
$(SIGNATURES)

Returns the estimated number of photoelectrons for a given ToT.
"""
function nphes(tot)
    if tot <= 20
        return 0.0
    end
    if tot <= 26
        return 1.0
    end
    if tot < 170
        return 1.0 + (tot - 26)/(1/0.28)
    end
    return 40.0 + (255 - tot)*2.0
end
