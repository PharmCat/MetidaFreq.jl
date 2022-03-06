

function freq(data, col; id = Dict())
    column = Tables.getcolumn(data, col)
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
