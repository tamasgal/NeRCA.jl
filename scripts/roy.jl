using NeRCA
using ROy
using Sockets

calib = Calibration("scripts/NeRCA_00000043_20022019.detx")
calib = Calibration("/home/tgal/data/detx/NeRCA_-00000001_20171212.detx")

client = NeRCA.CHClient(ip"127.0.0.1", 5553, ["IO_EVT"])

try
    for message in client
        event = read(IOBuffer(message.data), NeRCA.DAQEvent)
        hits = event.snapshot_hits
        chits = NeRCA.calibrate(calib, hits)
        thits = event.triggered_hits
        cthits = NeRCA.calibrate(calib, thits)
        println(length(cthits))
    end
catch InterruptException
    println("ok")
end
