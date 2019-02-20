using KM3NeT
using ROy
using Sockets

calib = KM3NeT.read_calibration("scripts/KM3NeT_00000043_20022019.detx")
calib = KM3NeT.read_calibration("/home/tgal/data/detx/KM3NeT_-00000001_20171212.detx")

client = KM3NeT.CHClient(ip"127.0.0.1", 5553, ["IO_EVT"])

try
    for message in client
        event = KM3NeT.read_io(IOBuffer(message.data), KM3NeT.DAQEvent)
        hits = event.snapshot_hits
        chits = KM3NeT.calibrate(hits, calib)
        thits = event.triggered_hits
        cthits = KM3NeT.calibrate(thits, calib)
        println(length(cthits))
    end
catch InterruptException
    println("ok")
end
