
# MetidaFreq.jl

struct ConTab <: AbstractData
    tab::Matrix{Int}
    rown::Union{Vector{String}, Nothing}
    coln::Union{Vector{String}, Nothing}
    id::Dict
end

"""
"""
function contab(m::AbstractMatrix{Int};
    rownames::Union{Vector{String}, Nothing} = nothing,
    colnames::Union{Vector{String}, Nothing} = nothing,
    id::Dict = Dict())
    ConTab(m, rown, coln, id)
end

"""
"""
function contab(data, row::Symbol, col::Symbol; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing, id = nothing)
    cols = Tables.columns(data)
    c    = isnothing(sort) ? (row, col) : isa(sort, Symbol) ? (row, col, sort) : Tuple(append!([row, col], sort))
    res  = contab_(Tuple(Tables.getcolumn(cols, y) for y in c))
    s    = size(res[1])
    if length(s) > 2
        dimn  = length(s) - 2
        ci    = CartesianIndices(s[3:end])
        v     = Vector{ConTab}(undef, length(ci))
        ckeys = Vector{Vector}(undef, dimn)
        for i = 3:length(s)
            ckeys[i] = collect(keys(res[2][i]))
        end
        for i = 1:length(ci)
            id = Dict{Symbol, promote_type(eltype.(ckeys)...)}()
            for j = 1:length(ci[i])
                id[c[2+j]] = ckeys[j][ci[i][j]]
            end
            v[i] = ConTab(Matrix(view(res[1], :, :, ci[i])), string.(collect(keys(res[2][1]))),  string.(collect(keys(res[2][2]))), id)
        end
        return DataSet(v)
    else
        return ConTab(res[1], string.(collect(keys(res[2][1]))),  string.(collect(keys(res[2][2]))), (isnothing(id) ? Dict() : id))
    end
end

function contab_(data::Tuple)
    d = Dict{Tuple{eltype.(data)...}, Int}()
    for (i, element) in enumerate(zip(data...))
        ind = ht_keyindex(d, element)
        if ind > 0
            @inbounds d.vals[ind] += one(Int)
        else
            @inbounds d[element] = one(Int)
        end
    end
    k    = collect(keys(d))
    n    = length(data)
    dims = Vector{Dict}(undef, n)
    for i = 1:n
        dims[i] = Dict{eltype(data[i]), Int}()
        s = Set{eltype(data[i])}()
        for j in 1:length(k)
            push!(s, k[j][i])
        end
        us = unique(s)
        for v = 1:length(us)
            dims[i][us[v]] = v
        end
    end
    m = zeros(Int, length.(dims)...)
    for j in 1:length(k)
        m[map(getindex, dims, k[j])...] = d[k[j]]
    end
    m, dims
end
