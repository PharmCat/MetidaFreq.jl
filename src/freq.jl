

function freq(data, col; id = Dict())
    if isa(col, String) cols = Symbol(col) else cols = col end
    column = Tables.getcolumn(data, cols)
    d = Dict{eltype(column), Int}()
    for i in column
        ind = ht_keyindex(d, i)
        if ind > 0
            @inbounds d.vals[ind] += one(Int)
        else
            @inbounds d[i] = one(Int)
        end
    end
    k = collect(keys(d))
    mx = Matrix{Int}(undef, 1, length(k))
    for i = 1:length(k)
        @inbounds mx[1, i] = d[k[i]]
    end
    ConTab(mx, [string(col)], string.(k), id)
end
