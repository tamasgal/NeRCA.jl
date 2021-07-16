"""
    function most_frequent(f::Function, iterable)

Return the most frequent value of a given iterable.
"""
function most_frequent(f::Function, iterable)
    d = Dict{Int, Int}()
    for element ∈ iterable
        v = f(element)
        if haskey(d, v)
            d[v] += 1
        else
            d[v] = 1
        end
    end
    candidate = 0
    key = 0
    for (k, v) in d
        if v > candidate
            key = k
            candidate = v
        end
    end
    return key
end


"""
    nthbitset(n, a) = !Bool((a >> (n - 1)) & 1)

Return `true` if the n-th bit of `a` is set, `false` otherwise.
"""
nthbitset(n, a) = Bool((a >> n) & 1)


"""
$(METHODLIST)

Categorise the struct elements of a vector by a given field into a dictionary of
`T.field => Vector{T}`.
"""
function categorize end

"""
$(TYPEDSIGNATURES)

Examples
========

```
julia> using NeRCA

julia> struct PMT
         dom_id
         time
       end

julia> pmts = [PMT(2, 10.4), PMT(4, 23.5), PMT(2, 42.0)];

julia> categorize(:dom_id, pmts)
Dict{Any, Vector{PMT}} with 2 entries:
  4 => [PMT(4, 23.5)]
  2 => [PMT(2, 10.4), PMT(2, 42.0)]
```
"""
@inline function categorize(field::Symbol, elements::Vector)
    _categorize(Val{field}(), elements)
end

"""
$(TYPEDSIGNATURES)
"""
@inline function _categorize(field::Val{F}, elements::Vector{T}) where {T,F}
    out = Dict{fieldtype(T, F), Vector{T}}()
    for el ∈ elements
        key = getfield(el, F)
        if !haskey(out, key)
            out[key] = T[]
        end
        push!(out[key], el)
    end
    out
end
