"""
    freqdict(data::AbstractVector) 

Return frequencies as `Dict`.
"""
function freqdict(data::AbstractVector) 
    d = Dict{eltype(data), Int}()
    for i in data
        ind = ht_keyindex(d, i)
        if ind > 0
            @inbounds d.vals[ind] += one(Int)
        else
            @inbounds d[i] = one(Int)
        end
    end
    d
end

"""
    freq(data, col; id = Dict())
 
Return frequencies as ix1 contingency table.
"""
function freq(data, col; id = Dict())
    if isa(col, String) cols = Symbol(col) else cols = col end
    column = Tables.getcolumn(data, cols)
    d = freqdict(column) 
    k = collect(keys(d))
    mx = Matrix{Int}(undef, length(k), 1 )
    for i = 1:length(k)
        @inbounds mx[i, 1] = d[k[i]]
    end
    ConTab(mx,  string.(k), [string(col)], id)
end

"""
    freq(data::AbstractVector; id = Dict()) 
 
Return frequencies as ix1 contingency table.
"""
function freq(data::AbstractVector; id = Dict()) 
    d = freqdict(data) 
    k = collect(keys(d))
    mx = Matrix{Int}(undef, length(k), 1)
    for i = 1:length(k)
        @inbounds mx[i, 1] = d[k[i]]
    end
    ConTab(mx, string.(k), [""], id)
end


