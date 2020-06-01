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
    event_hits = NeRCA.read_hits(NeRCA.OnlineFile(ROOTFILE))
    println(length(event_hits))
end

main()
