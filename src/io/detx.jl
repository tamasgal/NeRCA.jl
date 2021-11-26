function Calibration(detid::Integer, runid::Integer)
    ascii_data = detx(detid, runid)
    Calibration(IOBuffer(ascii_data))
end


function Calibration(filename::AbstractString)
    open(filename, "r") do fobj
        Calibration(fobj)
    end
end


function Calibration(io::IO)
    lines = readlines(io)
    filter!(e->!startswith(e, "#") && !isempty(strip(e)), lines)

    first_line = lowercase(first(lines))  # version can be v or V, halleluja

    if occursin("v", first_line)
        det_id, version = map(x->parse(Int,x), split(first_line, 'v'))
        n_modules = parse(Int, lines[4])
        idx = 5
    else
        det_id, n_modules = map(x->parse(Int,x), split(first_line))
        version = 1
        idx = 2
    end

    pos = Dict{Int32,Vector{NeRCA.Position}}()
    dir = Dict{Int32,Vector{NeRCA.Direction}}()
    t0s = Dict{Int32,Vector{Float64}}()
    dus = Dict{Int32,UInt8}()
    floors = Dict{Int32,UInt8}()
    omkeys = Dict{OMKey,Int32}()
    pmts = Dict{Int32,Tuple{Int32, UInt8}}()  # DOM ID and channel ID

    max_z = 0.0
    for mod ∈ 1:n_modules
        elements = split(lines[idx])
        module_id, du, floor = map(x->parse(Int, x), elements[1:3])
        if version == 4
            xₘ, yₘ, zₘ, q0, qx, qy, qz, t0ₘ = map(x->parse(Float64, x), elements[4:end-1])
        end
        n_pmts = parse(Int, elements[end])
        pos[module_id] = Vector{NeRCA.Position}()
        dir[module_id] = Vector{NeRCA.Direction}()
        t0s[module_id] = Vector{Float64}()
        dus[module_id] = du
        floors[module_id] = floor

        for pmt in 1:n_pmts
            l = split(lines[idx+pmt])
            pmt_id = parse(Int,first(l))
            x, y, z, dx, dy, dz = map(x->parse(Float64, x), l[2:7])
            max_z = max(max_z, z)
            t0 = parse(Float64,l[8])
            push!(pos[module_id], Position(x, y, z))
            push!(dir[module_id], Direction(dx, dy, dz))
            push!(t0s[module_id], t0)
            omkeys[OMKey(du, floor)] = module_id
            pmts[pmt_id] = (module_id, pmt-1)
        end
        idx += n_pmts + 1
    end
    n_dus = length(unique(values(dus)))
    Calibration(det_id, pos, dir, t0s, dus, floors, omkeys, pmts, max_z, n_dus)
end

@deprecate read_calibration(filename::AbstractString) Calibration(filename::AbstractString)
