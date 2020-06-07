println("Loading libraries...")
using NeRCA

if length(ARGS) < 2
    println("Usage: ./single_du_fit.jl DETX ROOTFILE")
    exit(1)
end

const ROOTFILE = ARGS[2]
# const calib = NeRCA.read_calibration("latest.detx")


function main()
    println("Starting reconstruction.")
    sparams = NeRCA.SingleDURecoParams()
    f = NeRCA.OnlineFile(ROOTFILE)
    event_shits = NeRCA.read_snapshot_hits(f)
    event_thits = NeRCA.read_triggered_hits(f)
    for (shits, thits) in zip(event_shits, event_thits)
        println(length(shits))
        println(length(thits))
        hits = NeRCA.combine(shits, thits)
        @show hits[1:4] length(hits)
    end
end

main()
