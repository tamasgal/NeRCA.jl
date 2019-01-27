function read_compound(dset::HDF5.HDF5Dataset, T::DataType)
    filetype = HDF5.datatype(dset) # packed layout on disk
    memtype_id = HDF5.h5t_get_native_type(filetype.id) # padded layout in memory
    @assert sizeof(T) == HDF5.h5t_get_size(memtype_id) "Type sizes don't match!"
    out = Vector{T}(undef, length(dset))
    HDF5.h5d_read(dset.id, memtype_id, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, out)
    HDF5.h5t_close(memtype_id)
    out
end


function read_compound(filename::AbstractString,
                       h5loc::AbstractString,
                       T::DataType)
    fobj = HDF5.h5open(filename)
    data = read_compound(fobj[h5loc], T)
    close(fobj)
    data
end


function read_hits(fobj::HDF5.HDF5File, idx::Int, n_hits::Int)
    hits = Vector{Hit}()
    channel_id = fobj["hits/channel_id"][idx+1:idx+n_hits]
    dom_id = fobj["hits/dom_id"][idx+1:idx+n_hits]
    t = fobj["hits/time"][idx+1:idx+n_hits]
    tot = fobj["hits/tot"][idx+1:idx+n_hits]
    triggered = fobj["hits/triggered"][idx+1:idx+n_hits]
    for i ∈ 1:n_hits
        hit =  Hit(channel_id[i], dom_id[i], t[i], tot[i], triggered[i])
        push!(hits, hit)
    end
    return hits
end


function read_hits(filename::AbstractString, event_id::Int)
    f = h5open(filename, "r")
    hit_indices = read_indices(f, "/hits")
    idx = hit_indices[event_id+1][1]
    n_hits = hit_indices[event_id+1][2]
    hits = read_hits(f, idx, n_hits)::Vector{Hit}
    close(f)
    return hits
end


function read_hits(filename::AbstractString,
                    event_ids::Union{Array{T}, UnitRange{T}}) where {T<:Integer}
    f = h5open(filename, "r")
    hit_indices = read_indices(f, "/hits")

    hits_collection = Dict{Int, Vector{Hit}}()
    for event_id ∈ event_ids
        idx = hit_indices[event_id+1][1]
        n_hits = hit_indices[event_id+1][2]
        hits = read_hits(f, idx, n_hits)::Vector{Hit}
        hits_collection[event_id] = hits
    end
    close(f)
    return hits_collection
end


function read_mchits(fobj::HDF5.HDF5File, idx::Int, n_hits::Int)
    hits = Vector{McHit}()
    a = fobj["mc_hits/a"][idx+1:idx+n_hits]
    origin = fobj["mc_hits/origin"][idx+1:idx+n_hits]
    pmt_id = fobj["mc_hits/pmt_id"][idx+1:idx+n_hits]
    t = fobj["mc_hits/time"][idx+1:idx+n_hits]
    for i ∈ 1:n_hits
        hit =  McHit(a[i], origin[i], pmt_id[i], t[i])
        push!(hits, hit)
    end
    return hits
end


function read_mchits(filename::AbstractString, event_id::Int)
    f = h5open(filename, "r")
    hit_indices = read_indices(f, "/mc_hits")
    idx = hit_indices[event_id+1][1]
    n_hits = hit_indices[event_id+1][2]
    hits = read_mchits(f, idx, n_hits)::Vector{McHit}
    close(f)
    return hits
end


function read_mchits(filename::AbstractString,
                    event_ids::Union{Array{T}, UnitRange{T}}) where {T<:Integer}
    f = h5open(filename, "r")
    hit_indices = read_indices(f, "/mc_hits")

    hits_collection = Dict{Int, Vector{McHit}}()
    for event_id ∈ event_ids
        idx = hit_indices[event_id+1][1]
        n_hits = hit_indices[event_id+1][2]
        hits = read_hits(f, idx, n_hits)::Vector{McHit}
        hits_collection[event_id] = hits
    end
    close(f)
    return hits_collection
end


function read_indices(filename::AbstractString, from::AbstractString)
    f = h5open(filename, "r")
    indices = read_indices(f, from)
    close(f)
    return indices
end


function read_indices(fobj::HDF5.HDF5File, from::AbstractString)
    idc = read(fobj, from * "/_indices")
    indices = [i.data for i ∈ idc]::Vector{Tuple{Int64,Int64}}
    return indices
end


function read_calibration(filename::AbstractString)
    lines = readlines(filename)

    if 'v' ∈ first(lines)
        det_id, version = map(x->parse(Int,x), split(first(lines), 'v'))
        n_doms = parse(Int, lines[4])
        idx = 5
    else
        det_id, n_doms = map(x->parse(Int,x), split(first(lines)))
        version = 1
        idx = 2
    end

    pos = Dict{Int32,Vector{KM3NeT.Position}}()
    dir = Dict{Int32,Vector{KM3NeT.Direction}}()
    t0s = Dict{Int32,Vector{Int32}}()
    dus = Dict{Int32,UInt8}()
    floors = Dict{Int32,UInt8}()

    for dom ∈ 1:n_doms
        dom_id, du, floor, n_pmts = map(x->parse(Int,x), split(lines[idx]))
        pos[dom_id] = Vector{KM3NeT.Position}()
        dir[dom_id] = Vector{KM3NeT.Direction}()
        t0s[dom_id] = Vector{Int32}()
        dus[dom_id] = du
        floors[dom_id] = floor

        for pmt in 1:n_pmts
            l = split(lines[idx+pmt])
            pmt_id = parse(Int,first(l))
            x, y, z, dx, dy, dz = map(x->parse(Float64, x), l[2:7])
            t0 = parse(Int,first(l[8]))
            push!(pos[dom_id], Position(x, y, z))
            push!(dir[dom_id], Direction(dx, dy, dz))
            push!(t0s[dom_id], t0)
        end
        idx += n_pmts + 1
    end

    Calibration(det_id, pos, dir, t0s, dus, floors)
end


function read_event_info(fobj::HDF5.HDF5File)
    event_info = Dict{Int32,EventInfo}()
    entries = read(fobj, "event_info")
    for entry in entries
        e = EventInfo(entry.data...)
        event_info[e.event_id] = e
    end
    event_info
end


function read_event_info(filename::AbstractString)
    event_info = h5open(filename) do file
        read_event_info(file)
    end
    event_info
end

function read_timeslice_info(fobj::HDF5.HDF5File)
    timeslice_info = Dict{Int32,TimesliceInfo}()
    entries = read(fobj, "event_info")
    for entry in entries
        e = TimesliceInfo(entry.data...)
        timeslice_info[e.event_id] = e
    end
    timeslice_info
end


function read_timeslice_info(filename::AbstractString)
    timeslice_info = h5open(filename) do file
        read_timeslice_info(file)
    end
    timeslice_info
end


#= struct EventReader =#
#=     filename::AbstractString =#
#=     n_events::Unsigned =#
#=     event_info::Dict{Int32,KM3NeT.EventInfo} =#
#=     hit_indices::Vector{Tuple{Int64,Int64}} =#
#=     tracks::Dict{Int32,Vector{MCTrack}} =#
#=     load_tracks::Bool =#
#=  =#
#=     EventReader(filename; load_tracks=false) = begin =#
#=         fobj = h5open(filename) =#
#=         n = length(read_event_info(fobj)) =#
#=         event_info = read_event_info(fobj) =#
#=         hit_indices = read_indices(fobj, "/hits") =#
#= #       mc_hit_indices = read_indices(fobj, "/mc_hits") =#
#=         if load_tracks =#
#=             tracks = read_tracks(fobj) =#
#=         else =#
#=             tracks = Dict{Int32,Vector{MCTrack}}() =#
#=         end =#
#=         close(fobj) =#
#=         new(filename, n, event_info, hit_indices, tracks, load_tracks) =#
#=     end =#
#= end =#


#= Base.start(E::EventReader) = begin =#
#=     fobj = h5open(E.filename) =#
#=     event_id = 0 =#
#=     return EventReaderState(fobj, event_id) =#
#= end =#
#=  =#
#=  =#
#= Base.next(E::EventReader, s) = begin =#
#=     event_id = s.event_id =#
#=     idx = E.hit_indices[s.event_id+1][1] =#
#=     n_hits = E.hit_indices[s.event_id+1][2] =#
#=     s.event_id += 1 =#
#=     if E.load_tracks =#
#=         tracks = E.tracks[event_id] =#
#=     else =#
#=         tracks = Vector{Track}() =#
#=     end =#
#=     event = Event(event_id, E.event_info[event_id], =#
#=                   read_hits(s.fobj, idx, n_hits), tracks) =#
#=     (event, s) =#
#= end =#
