@StructIO.io struct DAQPreamble
    length::Int32
    datatype::Int32
end

@StructIO.io struct DAQHeader
    preamble::DAQPreamble
    det_id::Int32
    run_id::Int32
    timeslice_id::Int32
    utc_seconds::UInt32
    ns_ticks::UInt32
end

@StructIO.io struct DAQEvent
    header::DAQHeader
    trigger_counter::UInt64
    trigger_mask::UInt64
    overlays::UInt32
end

@KM3NeT.StructIO.io struct TriggeredHit
    dom_id::Int32
    channel_id::Int8
    time::Int32
    tot::UInt8
    triggermask::UInt64
end
