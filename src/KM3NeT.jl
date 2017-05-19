module KM3NeT

using HDF5

export TimesliceHit, rows


abstract Hit

immutable RawHit <: Hit
    channel_id::UInt8
    time::Int32
    tot::UInt8
    dom_id::UInt32
end

RawHit(hit::HDF5.HDF5Compound{15}) = begin
    KM3NeT.RawHit(hit.data[1], hit.data[12], hit.data[13], hit.data[5])
end

immutable TimesliceHit <: Hit
    channel_id::Int8
    time::Int32
    tot::Int16
end


Base.show(io::IO, h::RawHit) = begin
    print(io, "$(typeof(h)): channel_id($(h.channel_id)), time($(h.time)), " *
              "tot($(h.tot)), dom_id($(h.dom_id))")
end


rows(x) = (x[i, :] for i in indices(x,1))


end # module
