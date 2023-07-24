Base.@kwdef struct MuonScanfitParameters
    tmax::Float64 = 25.0  # [ns]
    roadwidth::Float64 = 200.0  # [m]
    tmaxextra::Float64 = 18.0  # [ns]
    nmaxhits::Int = 50  # maximum number of hits to use
end

struct MuonScanfit
    params::MuonScanfitParameters
    detector::Detector
    directions::Vector{Direction{Float64}}
end
function Base.show(io::IO, m::MuonScanfit)
   print(io, "$(typeof(m)) in $(length(m.directions)) directions")
end

function (msf::MuonScanfit)(hits::Vector{T}) where T<:KM3io.AbstractHit
    l1builder = L1Builder(L1BuilderParameters(msf.params.tmax, false))
    rhits = l1builder(HitR1, msf.detector, hits)

    sort!(rhits)
    unique!(h->h.dom_id, rhits)

    clique = Clique(Match3B(msf.params.roadwidth, msf.params.tmaxextra))
    clusterize!(rhits, clique)



    clique1D = Clique(Match1D(msf.params.roadwidth, msf.params.tmaxextra))

    match_test = Match1D(msf.params.roadwidth, msf.params.tmax)


    for dir ∈ msf.directions
        rhits_copy = copy(rhits)
        R = rotator(dir)

        # rotate hits
        for (idx, rhit) ∈ enumerate(rhits_copy)
            rhits_copy[idx] = @set rhit.pos = R * rhit.pos
        end

        length(rhits_copy) > msf.params.nmaxhits && resize!(rhits_copy, msf.params.nmaxhits)
        sort!(rhits_copy; by=timetoz)
        println("...")
        @show length(rhits_copy)
        clusterize!(rhits_copy, clique1D)
        @show length(rhits_copy)

        length(rhits_copy) <= 3 continue  # TODO 3 comes from the number of parameters, retrieve from Line1Z fitter via type!

        # TODO x-y scane

    end
    rhits
end

struct Line1Z
    pos::Position{Float64}
    t::Float64
end
