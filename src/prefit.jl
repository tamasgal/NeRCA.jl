struct Hit
    t::Float64
    tot::Float64
end

starttime(hit) = hit.t
endtime(hit) = hit.t + hit.tot

function combine(hits::Vector{Hit})
    isempty(hits) && return Hit(0.0, 0.0)

    t1 = starttime(first(hits))
    t2 = endtime(first(hits))

    @inbounds for hit ∈ hits[2:end]
        _t1 = starttime(hit)
        _t2 = endtime(hit)
        if t1 > _t1
            t1 = _t1
        end
        if t2 < _t2
            t2 = _t2
        end
    end
    Hit(t1, t2 - t1)
end


