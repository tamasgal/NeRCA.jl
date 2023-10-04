"""
Performs the prefit algorithm which was used in DUMAND II.
"""
function dumandfit(hits::Vector{T}) where T <: AbstractCalibratedHit
    N = length(hits)
    D = 0.0
    dir = Direction(0.0, 0.0, 0.0)
    pos = Position(0.0, 0.0, 0.0)
    pes = [max(1, (h.tot - 26.3) / 4.5) for h in hits]
    # pes = [h.tot for h in hits]
    # pes = [h.multiplicity.count for h in hits]
    @inbounds for i ∈ 1:N
        @inbounds for k ∈ 1:N
            if i == k
                continue
            end
            t_i = hits[i].t
            t_k = hits[k].t
            q_ik = pes[i] * pes[k]
            t_ik = t_i * t_k
            pos += q_ik * (hits[i].pos*(t_k^2 - t_ik) + hits[k].pos*(t_i^2 - t_ik))
            dir += q_ik * (hits[i].pos - hits[k].pos) * (t_i - t_k)
            D += q_ik * (t_i - t_k)^2
        end
    end
    dir = normalize(dir/D)
    return Track(dir, pos/D, 0)
end
