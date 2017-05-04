module KM3NeT

export TimesliceHit, rows


abstract Hit

immutable TimesliceHit <: Hit
    channel_id::Int8
    time::Int32
    tot::Int16
end

rows(x) = (x[i, :] for i in indices(x,1))


end # module
