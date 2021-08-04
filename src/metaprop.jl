
struct MetaProp
    data::Vector{ConTab}
    metric::Symbol
    y::Vector
    var::Vector
end

struct MetaPropResult{Symbol}
    data::MetaProp
    wts::Vector{Float64}
    est::Float64
    var::Float64
    chisq::Float64
    hetq::Float64
    heti::Float64
    hettau::Float64
end

"""
    metaprop(d, metric; adj = 0)

Meta-analysis for 2x2 tables.

`ds`:

`metric`:
- :rr
- :or
- :diff
"""
function metaprop(d::DataSet, metric; adj = 0)
    if metric == :diff
        cty = contabdiff.(d.ds; adj = adj)
    elseif metric == :or
        cty = contabor.(d.ds; adj = adj)
    elseif metric == :rr
        cty = contabrr.(d.ds; adj = adj)
    else
    end
    y   = getindex.(cty, 1)
    var = getindex.(cty, 2)
    MetaProp(d.ds, metric, y, var)
end

"""
    metapropfixed(mp; weights = :default)

Inverce Variance method used by default.

`weights`:

- `:iv` | `:default`
- `:mh`

"""
function metapropfixed(mp; weights = :default)
    varwts    = 1 ./ mp.var
    if weights == :default || weights == :iv
        wts    = varwts
    elseif weights == :mh
        if mp.metric == :diff
            wts    = mhwdiff(mp.data)
        elseif mp.metric == :or
            wts    = mhwor(mp.data)
        elseif mp.metric == :rr
            wts    = mhwrr(mp.data)
        end
    #elseif weights == :peto
    else
        error("weights keyword unknown!")
    end
    k       = length(mp.y)
    var     = 1 / sum(varwts)
    est     = sum(wts .* mp.y) / sum(wts)
    chisq   = sum(varwts .* (mp.y .^ 2))
    q       = sum(varwts .* ((mp.y .- est) .^ 2))
    i²      = max(0, (q - (k - 1))/q * 100)
    MetaPropResult{:fixed}(mp, wts, est, var,  chisq, q, i², NaN)
end

"""
    metaproprandom(mp; tau = :default)

tau - τ² calculation method:
- `:dl` DerSimonian-Laird
- `:ho`
- `:hm` Hartung and Makambi
- `:sj` Sidik and Jonkman

"""
function metaproprandom(mp; tau = :default)
    k       = length(mp.y)
    wts     = 1 ./ mp.var
    swts    = sum(wts)
    var     = 1 / swts
    est     = sum(wts .* mp.y) / swts
    chisq   = sum(wts .* (mp.y .^ 2))
    q       = sum(wts .* ((mp.y .- est) .^ 2))
    i²      = max(0, (q - (k - 1))/q * 100)
    s       = sum(wts) - sum(wts .^ 2) / sum(wts)
    if tau == :dl || tau == :default # DerSimonian-Laird
        τ² = max(0, (q - (k - 1)) / s)
    elseif tau == :ho
        τ² = max(0, 1 / (k - 1) * sum((mp.y .- mean(mp.y)) .^ 2) - 1 / k * sum(mp.var))
    elseif tau == :hm # Hartung and Makambi
        τ² = q ^ 2 / (2 * (k - 1) + q) / (sum(wts) - sum(wts .^ 2) / sum(wts))
    elseif tau == :sj # Sidik and Jonkman
        qi  =  1 ./ (mp.var ./ max(0.01, 1 / (k - 1) * sum((mp.y .- mean(mp.y)) .^ 2) - 1 / k * sum(mp.var)) .+ 1)
        rest = sum(qi .* mp.y) / sum(qi)
        τ²   = 1 / (k - 1) * sum(qi .* ((mp.y .- rest) .^ 2))
    else
        error("tau keyword unknown!")
    end
    rwts = 1 ./ (mp.var .+ τ²)
    est     = sum(rwts .* mp.y) / sum(rwts)
    var     = 1 / sum(rwts) #?
    i²      = (τ² / (τ² + (k-1)/s)) * 100
    MetaPropResult{:random}(mp, rwts, est, var, chisq, q, i², τ²)
end

"""
    StatsBase.confint(mpr::MetaPropResult; level = 0.95)

Confidence interval for pooled proportion.
"""
function StatsBase.confint(mpr::MetaPropResult; level = 0.95)
    alpha = 1 - level
    se = sqrt(mpr.var)
    d = quantile(Normal(), 1 - alpha / 2) * se
    mpr.est - d, mpr.est + d
end

function mhwdiff(data)
    wts = Vector{Float64}(undef, length(data))
    for i = 1:length(data)
        wts[i] = (data[i].tab[1,1] + data[i].tab[1,2]) * (data[i].tab[2,1] + data[i].tab[2,2] ) / (data[i].tab[1,1] + data[i].tab[1,2] + data[i].tab[2,1] + data[i].tab[2,2])
    end
    wts
end
function mhwor(data)
    wts = Vector{Float64}(undef, length(data))
    for i = 1:length(data)
        wts[i] = data[i].tab[1,2] * data[i].tab[2,1] / (data[i].tab[1,1] + data[i].tab[1,2] + data[i].tab[2,1] + data[i].tab[2,2])
    end
    wts
end
function mhwrr(data)
    wts = Vector{Float64}(undef, length(data))
    for i = 1:length(data)
        wts[i] = (data[i].tab[1,1] + data[i].tab[1,2]) * data[i].tab[2,1] / (data[i].tab[1,1] + data[i].tab[1,2] + data[i].tab[2,1] + data[i].tab[2,2])
    end
    wts
end

#function petow(data)
#end
function contabdiff(contab; adj = 0)
    a = contab.tab[1,1] + adj
    b = contab.tab[1,2] + adj
    c = contab.tab[2,1] + adj
    d = contab.tab[2,2] + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return  pt - pc, pt * (1 - pt) / (a + b) + pc * (1 - pc) / (c + d)
end
function contabor(contab; adj = 0)
    a = contab.tab[1,1] + adj
    b = contab.tab[1,2] + adj
    c = contab.tab[2,1] + adj
    d = contab.tab[2,2] + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return log((pt / (1 - pt)) / (pc / (1 - pc))), 1 / a + 1 / b + 1 / c + 1 / d
end
function contabrr(contab; adj = 0)
    a = contab.tab[1,1] + adj
    b = contab.tab[1,2] + adj
    c = contab.tab[2,1] + adj
    d = contab.tab[2,2] + adj
    pt = a / (a + b)
    pc = c / (c + d)
    return log(pt / pc), 1 / a - 1 / (a + b) + 1 / c - 1 /(c + d)
end
