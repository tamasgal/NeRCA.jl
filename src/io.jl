struct OfflineFile end

struct OnlineFile
    fobj::UnROOT.ROOTFile

    OnlineFile(filename::AbstractString) = new(UnROOT.ROOTFile(filename))
end

struct KM3NETDAQSnapshotHit <: UnROOT.CustomROOTStruct
    dom_id::Int32
    channel_id::UInt8
    time::Int32
    tot::UInt8
end
function UnROOT.readtype(io, T::Type{KM3NETDAQSnapshotHit})
    T(UnROOT.readtype(io, Int32), read(io, UInt8), read(io, Int32), read(io, UInt8))
end

struct KM3NETDAQTriggeredHit <: UnROOT.CustomROOTStruct
    dom_id::Int32
    channel_id::UInt8
    time::Int32
    tot::UInt8
    trigger_mask::UInt64
end
UnROOT.packedsizeof(::Type{KM3NETDAQTriggeredHit}) = 24  # incl. cnt and vers

function UnROOT.readtype(io, T::Type{KM3NETDAQTriggeredHit})
    dom_id = UnROOT.readtype(io, Int32)
    channel_id = read(io, UInt8)
    tdc = read(io, Int32)
    tot = read(io, UInt8)
    cnt = read(io, UInt32)
    vers = read(io, UInt16)
    trigger_mask = UnROOT.readtype(io, UInt64)
    T(dom_id, channel_id, tdc, tot, trigger_mask)
end

function read_snapshot_hits(f::OnlineFile)
    data, offsets = UnROOT.array(f.fobj, "KM3NET_EVENT/KM3NET_EVENT/snapshotHits"; raw=true)
    UnROOT.splitup(data, offsets, KM3NETDAQSnapshotHit, skipbytes=10)
end

function read_triggered_hits(f::OnlineFile)
    data, offsets = UnROOT.array(f.fobj, "KM3NET_EVENT/KM3NET_EVENT/triggeredHits"; raw=true)
    UnROOT.splitup(data, offsets, KM3NETDAQTriggeredHit, skipbytes=10)
end

function Calibration(filename::AbstractString)
    lines = readlines(filename)
    filter!(e->!startswith(e, "#") && !isempty(strip(e)), lines)

    if 'v' ∈ first(lines)
        det_id, version = map(x->parse(Int,x), split(first(lines), 'v'))
        n_doms = parse(Int, lines[4])
        idx = 5
    else
        det_id, n_doms = map(x->parse(Int,x), split(first(lines)))
        version = 1
        idx = 2
    end

    pos = Dict{Int32,Vector{NeRCA.Position}}()
    dir = Dict{Int32,Vector{NeRCA.Direction}}()
    t0s = Dict{Int32,Vector{Float64}}()
    dus = Dict{Int32,UInt8}()
    floors = Dict{Int32,UInt8}()
    omkeys = Dict{Int32,OMKey}()

    max_z = 0.0
    for dom ∈ 1:n_doms
        dom_id, du, floor, n_pmts = map(x->parse(Int,x), split(lines[idx]))
        pos[dom_id] = Vector{NeRCA.Position}()
        dir[dom_id] = Vector{NeRCA.Direction}()
        t0s[dom_id] = Vector{Float64}()
        dus[dom_id] = du
        floors[dom_id] = floor

        for pmt in 1:n_pmts
            l = split(lines[idx+pmt])
            pmt_id = parse(Int,first(l))
            x, y, z, dx, dy, dz = map(x->parse(Float64, x), l[2:7])
            max_z = max(max_z, z)
            t0 = parse(Float64,l[8])
            push!(pos[dom_id], Position(x, y, z))
            push!(dir[dom_id], Direction(dx, dy, dz))
            push!(t0s[dom_id], t0)
            omkeys[pmt_id] = OMKey(dom_id, pmt-1)
        end
        idx += n_pmts + 1
    end
    n_dus = length(unique(values(dus)))
    Calibration(det_id, pos, dir, t0s, dus, floors, omkeys, max_z, n_dus)
end

@deprecate read_calibration(filename::AbstractString) Calibration(filename::AbstractString)

# Triggers
is_3dmuon(e::DAQEvent) = Bool(e.trigger_mask & 16 > 0)
is_3dshower(e::DAQEvent) = Bool(e.trigger_mask & 2 > 0)
is_mxshower(e::DAQEvent) = Bool(e.trigger_mask & 4 > 0)


# Online DAQ readout

function Base.read(s::IO, ::Type{T}) where T<:DAQEvent
    length = read(s, Int32)
    type = read(s, Int32)
    det_id = read(s, Int32)
    run_id = read(s, Int32)
    timeslice_id = read(s, Int32)
    timestamp = read(s, Int32)
    ticks = read(s, Int32)
    trigger_counter = read(s, Int64)
    trigger_mask = read(s, Int64)
    overlays = read(s, Int32)

    n_triggered_hits = read(s, Int32)
    triggered_hits = Vector{TriggeredHit}()
    sizehint!(triggered_hits, n_triggered_hits)
    triggered_map = Dict{Tuple{Int32, UInt8, Int32, UInt8}, Int64}()
    @inbounds for i ∈ 1:n_triggered_hits
        dom_id = read(s, Int32)
        channel_id = read(s, UInt8)
        time = bswap(read(s, Int32))
        tot = read(s, UInt8)
        trigger_mask = read(s, Int64)
        triggered_map[(dom_id, channel_id, time, tot)] = trigger_mask
        push!(triggered_hits, TriggeredHit(dom_id, channel_id, time, tot, trigger_mask))
    end

    n_hits = read(s, Int32)
    hits = Vector{Hit}()
    sizehint!(hits, n_hits)
    @inbounds for i ∈ 1:n_hits
        dom_id = read(s, Int32)
        channel_id = read(s, UInt8)
        time = bswap(read(s, Int32))
        tot = read(s, UInt8)
        key = (dom_id, channel_id, time, tot)
        triggered = false
        if haskey(triggered_map, key)
            triggered = true
        end
        push!(hits, Hit(channel_id, dom_id, time, tot, triggered))
    end

    T(det_id, run_id, timeslice_id, timestamp, ticks, trigger_counter, trigger_mask, overlays, n_triggered_hits, triggered_hits, n_hits, hits)
end
