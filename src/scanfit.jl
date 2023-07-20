struct MuonScanfitParameters
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    directions::Vector{Direction{Float64}}
end
function Base.show(io::IO, mp::MuonScanfit)
   print(io, "$(typeof(mp)) in $(length(mp.directions)) directions")
end

function (mp::MuonScanfit)(hits)
    # TODO: make it user-definable
    l1builder = L1Builder(L1BuilderParameters(25, false))
    rhits = l1builder(HitR1, mp.detector, hits)
    sort!(rhits)
    unique!(h->h.dom_id, rhits)

    for dir ∈ mp.directions
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits)
            rhits[idx] = @set rhit.pos = R * rhit.pos
        end



        # rotate hits back
        for (idx, rhit) ∈ enumerate(rhits)
            rhits[idx] = @set rhit.pos = R \ rhit.pos
        end
    end
    rhits
end

struct Line1Z
    pos::Position{Float64}
    t::Float64
end
