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
    rhits = calibrate(HitR1, mp.detector, hits)

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

struct L2BuilderParameters
    n_hits::Int
    Δt::Float64
    ctmin::Float64
end

struct L2Builder
    params::L2BuilderParameters
end
