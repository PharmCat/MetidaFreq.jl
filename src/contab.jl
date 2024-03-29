
# MetidaFreq.jl

struct ConTab{T <: AbstractMatrix{Int}, R <: Vector{S} where S <: AbstractString, C <:  Vector{S} where S <: AbstractString, ID <: Dict} <: AbstractIdData
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
    if isnothing(rownames)
        rownames = fill("", size(m, 1))
    end
    if isnothing(colnames)
        colnames = fill("", size(m, 2))
    end
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
    if isnothing(rownames)
        rownames = [""]
    end
    if isnothing(colnames)
        colnames = fill("", length(v))
    end
    contab(permutedims(v); rownames = rownames, colnames = colnames, id = id)
end

"""
    contab(ct::ConTab, rr, cr)

Make ConTab with `ct`, rows `rr` and columns `cr`.
"""
function contab(ct::ConTab, rr, cr)
    if isa(rr, Int) rr = [rr] end
    if isa(cr, Int) cr = [cr] end
    contab(ct.tab[rr, cr]; rownames = ct.rown[rr], colnames = ct.coln[cr], id = copy(ct.id))
end

"""
    Base.getindex(ct::ConTab, args...) = getindex(ct.tab, args...)

Returns values of Contab by index.
"""
Base.getindex(ct::ConTab, args...) = getindex(ct.tab, args...)

"""
    Base.permutedims(ct::ConTab)

ConTab permutedims.
"""
function Base.permutedims(ct::ConTab)
    ConTab(permutedims(ct.tab), copy(ct.coln), copy(ct.rown), copy(ct.id))
end

function sumrows_(f::Function, ct::ConTab)
    n  = size(ct.tab, 1)
    mx = Matrix{Int}(undef, n, 1)
    for i = 1:n
        mx[i,1] = sum(f, view(ct.tab, i, :))
    end
    mx
end
"""
    sumrows(f::Function, contab::ConTab; coln = "Val")

Aplpy function to each element of row, sum and make new column.
"""
function sumrows(f::Function, ct::ConTab; coln = "Val")
    mx = sumrows_(f, ct)
    contab(mx;
        rownames = copy(ct.rown),
        colnames = [coln],
        id       = copy(ct.id))
end
"""
    sumrows(contab::ConTab; coln = "Val")

Aplpy `identity` function to each element of row, sum and make new column.
"""
function sumrows(ct::ConTab; coln = "Val")
    sumrows(identity, ct; coln = coln)
end

"""
    addcol(ct::ConTab, col::Vector{Int}; coln = "Val")

Add column.
"""
function addcol(ct::ConTab, col::Vector{Int}; coln = "Val")
    contab(hcat(ct.tab, col);
        rownames = copy(ct.rown),
        colnames = push!(copy(ct.rown), coln),
        id       = copy(ct.id))
end

"""
    addcol(f::Function, ct::ConTab; coln = "Val")

Apply function to row and make new column.
"""
function addcol(f::Function, ct::ConTab; coln = "Val")
    n   = size(ct, 1)
    col = Vector{Int}(undef, n)
    for i = 1:n
        col[i] = f(view(ct.tab, i, :))
    end
    contab(hcat(ct.tab, col);
        rownames = copy(ct.rown),
        colnames = push!(copy(ct.rown), coln),
        id       = copy(ct.id))
end

"""
    addcol(f::Function, ct::ConTab, col::Vector{Int}; coln = "Val")

Example function (x,y) -> sum(x) + y, where x - row, y - value of col item. 
"""
function addcol(f::Function, ct::ConTab, col::Vector{Int}; coln = "Val")
    n   = size(ct, 1)
    ncol = Vector{Int}(undef, n)
    for i = 1:n
        ncol[i] = f(view(ct.tab, i, :), col[i])
    end
    contab(hcat(ct.tab, ncol);
        rownames = copy(ct.rown),
        colnames = push!(copy(ct.coln), coln),
        id       = copy(ct.id))
end
"""
    colreduce(f::Function, data::DataSet{<:ConTab}; coln = nothing)

Sum rows for each table, than make new table where in each column pleced sums.
"""
function colreduce(f::Function, data::DataSet{<:ConTab}; coln = nothing)
    fst = data.ds[1].rown
    if length(data) > 1
        for i = 2:length(data)
            if fst != data.ds[i].rown
                error("row names not equal")
            end
        end
    end
    if isnothing(coln)
        coln = Vector{String}(undef, length(data))
        for i = 1:length(data)
            coln[i] = string(data.ds[i].id)
        end
    end
    mx = hcat([sumrows_(f, i) for i in data.ds]...)
    contab(mx;
        rownames = copy(fst),
        colnames = coln,
        id       = Dict())
end

"""
    colorder(ct::ConTab, v::AbstractVector{AbstrctString})

Make contigency trable with column order as in `v` vector. If item in `v` not present in tab collumn's names - collumn of zeros will be made for that item.

Example 

```
Contingency table:
--------- -------- ----- -------
              B      D    Total 
--------- -------- ----- -------
   G1         1     123     124
   G2         0     124     124
--------- -------- ----- -------
 ID: ColName => id;
```

after `colorder(ct, ["A", "B", "C", "D"])`:

```
Contingency table:
--------- -------- --------- --------- ----- -------
              A         B         C      D    Total 
--------- -------- --------- --------- ----- -------
    G1        0         1         0     123     124
    G2        0         0         0     124     124
--------- -------- --------- --------- ----- -------
ID: ColName => id;
```
"""
function colorder(ct::ConTab, v::AbstractVector{<:AbstractString}) 

    foreach(x-> if x ∉ v error("Not all col names in new vector: "*x*" ID: "*string(ct.id)) end, ct.coln)
     
    tab = zeros(Int, size(ct.tab, 1), length(v))
    for i = 1:length(v)
        ind = findfirst(x -> x == v[i], ct.coln)
        if !isnothing(ind) tab[:, i] .= view(ct.tab, :, ind) end
    end
    ConTab(tab, ct.rown, v, ct.id)
end

function colorder(data::DataSet{<:ConTab}, v::AbstractVector{<:AbstractString})
    DataSet(map(x -> colorder(x, v), getdata(data)))
end

function contab(data, row::Union{Symbol, String}, col; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing)
    namevec = names(data)
    if eltype(col) <: Int icol = namevec[first(col)] else icol = first(col) end
    tab = contab(data, row, icol; sort = sort, idp = :ColName => icol)
    if isa(tab, ConTab) ds = DataSet([tab]) else ds = tab end
    if length(col) > 1
        for i = 2:length(col)
            if isa(col[i], Int) icol = namevec[col[i]] else icol = col[i] end
            tab = contab(data, row, icol; sort = sort, idp = :ColName => icol)
            if isa(tab, ConTab) push!(getdata(ds), tab) else append!(getdata(ds), getdata(tab))  end
        end
    end
    return ds
end
"""
    contab(data, row::Symbol, col::Symbol; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing, id = nothing)

Make contingency table from data using `row` and `col` columns.
"""
function contab(data, row::Union{Symbol, String}, col::Union{Symbol, String}; sort::Union{Nothing, Symbol, AbstractVector{Symbol}} = nothing, idp::Union{Nothing, Pair} = nothing)
    if isa(row, String) row = Symbol(row) end
    if isa(col, String) col = Symbol(col) end
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
                rowstr[rdi[k]] = string(k)
            end
            for k in keys(cdi)
                colstr[cdi[k]] = string(k)
            end
            if !isnothing(idp) push!(id, idp) end
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
        return ConTab(res[1], rowstr,  colstr, (isnothing(idp) ? Dict() : Dict(idp)))
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

"""
    dropzeros!(ds::DataSet{<:ConTab})

Drop tables from dataset if no observation.
"""
function dropzeros!(ds::DataSet{<:ConTab})
    inds = Int[]
    for i in 1:length(ds)
        if sum(ds[i].tab) == 0 push!(inds, i) end
    end
    if length(inds) > 0
        deleteat!(ds.ds, inds)
    end
    ds
end

function Base.size(contab::ConTab)
    size(contab.tab)
end
function Base.size(contab::ConTab, dim::Int)
    size(contab.tab, dim)
end
function Base.show(io::IO, contab::ConTab)
    println(io, "  Contingency table:")
    tab  = hcat(contab.tab, sum(contab.tab, dims = 2))
    coln = Vector{String}(undef, length(contab.coln) + 1)
    copyto!(view(coln, 1:length(contab.coln)), contab.coln)
    coln[end] = "Total"
    #coln = push!(copy(contab.coln), "Total")
    PrettyTables.pretty_table(io, tab; header = coln, row_labels = contab.rown, tf = PrettyTables.tf_compact)
    if !isnothing(contab.id) && length(contab.id) > 0
        print(io, "  ID: ")
        for (k,v) in contab.id
            print(io, "$k => $v; ")
        end
    end
end

function Base.show(io::IO, ds::DataSet{<:ConTab})
    for i in 1:length(ds)
        println(io, ds[i])
    end
end


function MetidaBase.metida_table_(obj::DataSet{T}) where T <: ConTab
    idset  = Set(keys(first(obj).id))
    s      = size(first(obj).tab)
    r      = first(obj).rown
    c      = first(obj).coln
    if length(obj) > 1
        for i = 2:length(obj)
            union!(idset,  Set(keys(obj[i].id)))
            if size(obj[i]) != s error("Unequal table size.") end
            if obj[i].rown != r error("Unequal row names.") end
            if obj[i].coln != c error("Unequal col names.") end
        end
    end
    # Check all same size
    
    inds = [(x,y) for x in 1:s[1] for y in 1:s[2]] 

    mt1 = metida_table_((getid(obj, :, c) for c in idset)...; names = idset)

    names = ["r$(i[1])("*first(obj).rown[i[1]]*"):c$(i[2])("*first(obj).coln[i[2]]*")" for i in inds]

    mt2 = metida_table_((map(x->x.tab[i[1], i[2]], getdata(obj)) for i in inds)...; names = names)
    merge(mt1, mt2)
end
