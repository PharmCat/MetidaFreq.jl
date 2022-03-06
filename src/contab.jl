
# MetidaFreq.jl

struct ConTab{T <: AbstractMatrix{Int}, R <: Union{Nothing, Vector{S}} where S <: AbstractString, C <: Union{Nothing, Vector{S}} where S <: AbstractString, ID <: Dict} <: AbstractIdData
    tab::T
    rown::R
    coln::C
    id::ID
end

"""
    contab(m::AbstractMatrix{Int};
        rownames::Union{Vector{String}, Nothing} = nothing,
        colnames::Union{Vector{String}, Nothing} = nothing,
        id::Dict = Dict())
"""
function contab(m::AbstractMatrix{Int};
    rownames::Union{Vector{String}, Nothing} = nothing,
    colnames::Union{Vector{String}, Nothing} = nothing,
    id::Dict = Dict())
    ConTab{typeof(m), typeof(rownames), typeof(colnames), typeof(id)}(m, rownames, colnames, id)
end
"""
    contab(v::AbstractVector{Int};
        rownames::Union{Vector{String}, Nothing} = nothing,
        colnames::Union{Vector{String}, Nothing} = nothing,
        id::Dict = Dict())
"""
function contab(v::AbstractVector{Int};
    rownames::Union{Vector{String}, Nothing} = nothing,
    colnames::Union{Vector{String}, Nothing} = nothing,
    id::Dict = Dict())
    contab(permutedims(v); rownames = rownames, colnames = colnames, id = id)
end

"""
    contab(ct::ConTab, rr, cr)

Make ConTab with `ct`, rows `rr` and columns `cr`.
"""
function contab(ct::ConTab, rr, cr)
    if isa(rr, Int) rr = [rr] end
    if isa(cr, Int) cr = [cr] end
    contab(ct.tab[rr, cr]; rownames = ct.rown[rr], colnames = ct.coln[cr], id = ct.id)
end

"""
    Base.permutedims(ct::ConTab)

ConTab permutedims.
"""
function Base.permutedims(ct::ConTab)
    ConTab(permutedims(ct.tab), ct.coln, ct.rown, ct.id)
end


"""
    contab(data, row::Symbol, col::Symbol; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing, id = nothing)
"""
function contab(data, row::Symbol, col::Symbol; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing, id = nothing)
    cols = Tables.columns(data)
    c    = isnothing(sort) ? (row, col) : isa(sort, Symbol) ? (row, col, sort) : Tuple(append!([row, col], sort))
    res  = contab_(Tuple(Tables.getcolumn(cols, y) for y in c))
    s    = size(res[1])
    rdi = res[2][1]
    cdi = res[2][2]
    rowstr = Vector{String}(undef, length(rdi))
    colstr = Vector{String}(undef, length(cdi))
    if length(s) > 2
        dimn  = length(s) - 2
        ci    = CartesianIndices(s[3:end]) #?
        v     = Vector{ConTab}(undef, length(ci))
        rdict = Vector{Dict}(undef, dimn)
        for i = 1:length(s)-2
            rdict[i] = Dict(v => k for (k,v) in res[2][2+i])
        end
        for i = 1:length(ci)
            id = Vector{Pair}(undef, length(ci[i]))
            for j = 1:length(ci[i])
                id[j] = c[2+j] => rdict[j][ci[i][j]]
            end
            for k in keys(rdi)
                rowstr[rdi[k]] = String(k)
            end
            for k in keys(cdi)
                colstr[cdi[k]] = String(k)
            end
            v[i] = ConTab(Matrix(view(res[1], :, :, ci[i])), rowstr,  colstr, Dict(id))
        end
        return DataSet(v)
    else
        for k in keys(rdi)
            rowstr[rdi[k]] = string(k)
        end
        for k in keys(cdi)
            colstr[cdi[k]] = string(k)
        end
        return ConTab(res[1], rowstr,  colstr, (isnothing(id) ? Dict() : id))
    end
end

function contab_(data::Tuple, T::Type = promote_type(eltype.(data)...))
    #d = Dict{Tuple{eltype.(data)...}, Int}()
    d = Dict{Tuple{map(eltype, data)...}, Int}()
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
    dims = Vector{Dict{T, Int}}(undef, n)
    @inbounds for i in Base.OneTo(n)
        dims[i] = Dict{eltype(data[i]), Int}()
        s = Set{eltype(data[i])}()
        @inbounds for j in Base.OneTo(length(k))
            push!(s, k[j][i])
        end
        us = collect(s)
        @inbounds for v in Base.OneTo(length(us))
            dims[i][us[v]] = v
        end
    end

    m = zeros(Int, length.(dims)...)
    @inbounds for j in Base.OneTo(length(k))
        m[map(getindex, dims, k[j])...] = d[k[j]]
    end

    #=
    diml = Tuple(map(length, dims))
    li   = LinearIndices(diml)
    m = zeros(Int, length(li))
    @inbounds for j in Base.OneTo(length(k))
        m[li[map(getindex, dims, k[j])...]] = d[k[j]]
    end
    reshape(m, diml...), dims
    =#
    m, dims
end

function getinvindex(data)
    data.pool.invindex
end

function contab_(data::NTuple{n, AbstractCategoricalVector}) where n

    levs = map(levels, data)
    dims = Tuple(map(length, levs))

    #ci = CartesianIndices(dims)
    #m  = zeros(Int, dims...)
    li = LinearIndices(dims)
    a = zeros(Int, length(li))

    r = Vector{Int}(undef, n)
    @inbounds for i in 1:length(data[1])
        @inbounds for j = 1:n
            r[j] = Int(data[j].refs[i])
        end
        a[li[r...]] += 1
        #m[r...] += 1
    end
    #m = reshape(a, dims...)
    reshape(a, dims...), getinvindex.(data)
end

function Base.size(contab::ConTab)
    size(contab.tab)
end
function Base.size(contab::ConTab, dim::Int)
    size(contab.tab, dim)
end
function Base.show(io::IO, contab::ConTab)
    println(io, "  Contingency table:")
    PrettyTables.pretty_table(io, contab.tab; header = contab.coln, row_names = contab.rown, tf = PrettyTables.tf_compact)
    if !isnothing(contab.id) && length(contab.id) > 0
        print(io, "  ID: ")
        for (k,v) in contab.id
            print(io, "$k => $v; ")
        end
    end
end
