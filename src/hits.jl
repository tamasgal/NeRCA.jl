"""
    function triggered(hits::Vector{T}) where {T<:Hit}

Return a `Vector` of triggered hits.
"""
triggered(hits::Vector{T}) where {T<:DAQHit} = filter(h->h.triggered, hits)


"""
    function nfoldhits(hits::Vector{T}, Δt, n) where {T<:Hit}

Create a `Vector` with hits contributing to `n`-fold coincidences within a time
window of Δt.
"""
function nfoldhits(hits::Vector{T}, Δt, n) where {T<:DAQHit}
    hit_map = DefaultDict{Integer}{Vector{T}}(() -> T[])
    for hit ∈ sort(hits)
        push!(hit_map[hit.dom_id], hit)
    end
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
    function domhits(hits::Vector{T}) where {T<:Hit}

Sort hits by DOM ID and put them into a dictionary.
"""
function domhits(hits::Vector{T}) where {T<:DAQHit}
    hit_map = DefaultDict{Integer}{Vector{T}}(() -> T[])
    for hit ∈ hits
        push!(hit_map[hit.dom_id], hit)
    end
    hit_map
end


"""
    function duhits(hits::Vector{T}) where {T<:CalibratedHit}

Sort hits by DU and put them into a dictionary.
"""
function duhits(hits::Vector{T}) where {T<:CalibratedHit}
    hit_map = DefaultDict{Integer}{Vector{T}}(() -> T[])
    for hit ∈ hits
        push!(hit_map[hit.du], hit)
    end
    hit_map
end
