
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

Meta-analysis for 2x2 tables. Where:

`d`: `DataSet{ConTab}

`metric`:

- :rr (Risk Ratio)
- :or (Odd Ratio)
- :diff (Risk Difference)

`adj` - adjustment value.

"""
function metaprop(d::DataSet, metric; adj = 0)
    if metric == :diff
        cty = contabdiff.(d.ds; adj = adj)
    elseif metric == :or
        cty = contabor.(d.ds; adj = adj)
    elseif metric == :rr
        cty = contabrr.(d.ds; adj = adj)
    else
        error("metric unknown")
    end
    y   = getindex.(cty, 1)
    var = getindex.(cty, 2)
    MetaProp(d.ds, metric, y, var)
end

"""
    metapropfixed(mp; weights = :default)

Inverce Variance method used by default.

`weights`:

- `:iv` | `:default` (Inverce Variance)
- `:mh` (Mantel Haenszel)

For Risk Difference `Sato, Greenland, & Robins (1989)` modification for variance estimation used. 

*Results for RR and OR are in log-scale.**
"""
function metapropfixed(mp; weights = :default)
    varwts    = 1 ./ mp.var
    if weights == :default || weights == :iv
        wts    = varwts
        est     = sum(wts .* mp.y) / sum(wts)
        var    = 1 / sum(varwts)
    elseif weights == :mh
        if mp.metric == :diff
            wts      = mhwdiff(mp.data)
            est      = sum(wts .* mp.y) / sum(wts)
            var      = mhvardiff(est, mp.data)
        elseif mp.metric == :or
            wts      = mhwor(mp.data)
            #est    = sum(wts .* mp.y) / sum(wts)
            est, var = mhvarestor(mp.data) 
        elseif mp.metric == :rr
            wts      = mhwrr(mp.data)
            #est     = sum(wts .* mp.y) / sum(wts)
            est, var = mhvarestrr(mp.data)
        end
    else
        error("weights keyword unknown!")
    end
    k       = length(mp.y)
    chisq   = sum(varwts .* (mp.y .^ 2))
    q       = sum(varwts .* ((mp.y .- est) .^ 2))
    i²      = max(0, (q - (k - 1))/q * 100)
    MetaPropResult{:fixed}(mp, wts, est, var,  chisq, q, i², NaN)
end

sq(x) = x * x
sumsq(x) = sum(sq, x)
"""
    metaproprandom(mp; tau = :default)

tau - τ² calculation method:

- `:dl` DerSimonian-Laird (by default)
- `:ho` Hedges - Olkin
- `:hm` Hartung and Makambi (Veroniki et al. 2016)
- `:sj` Sidik and Jonkman
- `:ml` Maximum likelihood (ML) method (Veroniki et al. 2016)
- `:reml` Restricted maximum likelihood (REML) method (Veroniki et al. 2016)

"""
function metaproprandom(mp; tau = :default)
    k       = length(mp.y)
    wts     = 1 ./ mp.var
    swts    = sum(wts)
    var     = 1 / swts
    est     = dot(wts, mp.y) / swts
    zipwtsy = zip(wts, mp.y)
    chisq   = sum(x-> x[1] * x[2]^2, zipwtsy)
    #chisq   = sum(wts .* (mp.y .^ 2))
    q       = sum(x-> x[1] * (x[2] - est)^2, zipwtsy)
    #q       = sum(wts .* ((mp.y .- est) .^ 2))
    i²      = max(0, (q - (k - 1))/q * 100)
    s       = swts - sumsq(wts) / swts
    #s       = sum(wts) - sum(wts .^ 2) / sum(wts)
    if tau == :dl || tau == :default # DerSimonian-Laird
        τ² = max(0, (q - (k - 1)) / s)
    elseif tau == :ho # Hedges - Olkin
        hom = mean(mp.y)
        ho = sum(x-> (x - hom)^2, mp.y) 
        τ² = max(0, 1 / (k - 1) * ho - 1 / k * sum(mp.var))
        #τ² = max(0, 1 / (k - 1) * sum((mp.y .- mean(mp.y)) .^ 2) - 1 / k * sum(mp.var))
    elseif tau == :hm # Hartung and Makambi
        τ² = q ^ 2 / (2 * (k - 1) + q) / s
    elseif tau == :sj # Sidik and Jonkman
        qi  =  1 ./ (mp.var ./ max(0.01, 1 / (k - 1) * sum((mp.y .- mean(mp.y)) .^ 2) - 1 / k * sum(mp.var)) .+ 1)
        rest = sum(qi .* mp.y) / sum(qi)
        τ²   = 1 / (k - 1) * sum(qi .* ((mp.y .- rest) .^ 2))
    elseif tau == :ml
        τ² = max(0, (q - (k - 1)) / s)
        τ² = taumlsolve(τ², k, mp.y, mp.var)
    elseif tau == :reml
        τ² = max(0, (q - (k - 1)) / s)
        τ² = tauremlsolve(τ², k, mp.y, mp.var)
    else
        error("tau keyword unknown!")
    end
    rwts = @. 1 / (mp.var + τ²)
    #rwts = 1 ./ (mp.var .+ τ²)
    est     = dot(rwts, mp.y) / sum(rwts)
    #est     = sum(rwts .* mp.y) / sum(rwts)
    var     = 1 / sum(rwts) 
    i²      = (τ² / (τ² + (k - 1) / s)) * 100
    MetaPropResult{:random}(mp, rwts, est, var, chisq, q, i², τ²)
end

"""
    StatsBase.confint(mpr::MetaPropResult; level = 0.95)

Confidence interval for pooled proportion.

!!! Warn
    Results are in log-scale for OR and RR. 
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

# Greenland & Robins (1985)
#=
function mhvardiff(data)
    var = 0.0
    d   = 0.0
    for i = 1:length(data)
        N  = data[i].tab[1,1] + data[i].tab[1,2] + data[i].tab[2,1] + data[i].tab[2,2]
        n1 = data[i].tab[1,1] + data[i].tab[1,2]
        n2 = data[i].tab[2,1] + data[i].tab[2,2]

        var += data[i].tab[1,1] / N^2 * data[i].tab[1,2] * n2^2 / n1 
        +      data[i].tab[2,1] / N^2 * data[i].tab[2,2] * n1^2 / n2
        d   += n1 * n2 / N
    end
    var / d^2
end
=#
# Sato, Greenland, & Robins (1989)
function mhvardiff(est, data)
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    for i = 1:length(data)
        a = data[i].tab[1,1]
        b = data[i].tab[1,2]
        c = data[i].tab[2,1]
        d = data[i].tab[2,2]
        n1 = a + b
        n2 = c + d
        N  = n1 + n2
        sum1 += c * (n1 / N) ^ 2 - a * (n2 / N) ^2 + (n1 / N) * (n2 / N) * (n2 - n1) / 2
        sum2 += a * d / N + c * b / N
        sum3 += n1 * n2 / N
    end
    (est * sum1 + sum2 / 2) / sum3 .^ 2
end
# variance by Robins et al. (1986)
function mhvarestor(data)
    Rsum  = 0.0
    Ssum  = 0.0
    sum1  = 0.0
    sum2  = 0.0
    sum3  = 0.0
    for i = 1:length(data)
        a = data[i].tab[1,1]
        b = data[i].tab[1,2]
        c = data[i].tab[2,1]
        d = data[i].tab[2,2]
        n1 = a + b
        n2 = c + d
        N  = n1 + n2
        P  = a/N + d/N
        Q  = b/N + c/N
        R  = (a/N) * d
        S  = (b/N) * c
        Rsum += R
        Ssum += S
        sum1 += P * R
        sum2 += P * S + Q * R
        sum3 += Q * S
    end
    log(Rsum / Ssum), (sum1 / Rsum^2 + sum2 / (Rsum * Ssum)  + sum3 / Ssum^2) / 2
end

function mhvarestrr(data)
    Rsum  = 0.0
    Ssum  = 0.0
    sum1  = 0.0
    for i = 1:length(data)
        a = data[i].tab[1,1]
        b = data[i].tab[1,2]
        c = data[i].tab[2,1]
        d = data[i].tab[2,2]
        n1 = a + b
        n2 = c + d
        N  = n1 + n2
        Rsum += a * (n2 / N) 
        Ssum += c * (n1 / N) 
        sum1 += (n1 / N) * (n2 / N) * (a + c) - (a / N) * c
    end
    log(Rsum / Ssum), sum1 / (Rsum * Ssum)
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

function tauml(μ, τ², k, y::AbstractVector, v::AbstractVector)
    sum1 = sum(x-> log(x + τ²), v)
    sum2 = sum(x-> (x[1] - μ)^2 /(x[2] + τ²), zip(y, v))
    - (k * log(2π) + sum1 + sum2) / 2
end

function taureml(μ, τ², k, y::AbstractVector, v::AbstractVector)
    tauml(μ, τ², k, y, v) - log(sum(x-> 1/(x + τ²), v))/2
end

function taumlsolve(τ², k, y::AbstractVector, v::AbstractVector)
    varwts = @. 1/(v + τ²)
    μ      =  dot(varwts, y)/sum(varwts)
    optf(x) = - tauml(x[1], exp(x[2]), k, y, v)
    res = optimize(optf, [μ, log(τ²)], Newton(); autodiff = :forward)
    mz = Optim.minimizer(res)
    exp(mz[2])
end

function tauremlsolve(τ², k, y::AbstractVector, v::AbstractVector)
    τ²     = taumlsolve(τ², k, y, v)
    varwts = @. 1/(v + τ²)
    μ      =  dot(varwts, y)/sum(varwts)
    optf(x) = - taureml(μ, exp(x[1]), k, y, v)
    res = optimize(optf, [log(τ²)], Newton(); autodiff = :forward)
    mz = Optim.minimizer(res)
    exp(mz[1])
end


function Base.show(io::IO, mp::MetaProp)
    println(io, "  Meta-proportion:")
    println(io, "  Tables: $(length(mp.data))")
    println(io, "  Metric: $(mp.metric)")
    println(io, "  Metric vector: $(mp.y)")
    print(io,   "  Metric variance: $(mp.var)")
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


function weights(mpr::MetaPropResult)
    mpr.wts ./ (sum(mpr.wts) / 100)
end

function Base.show(io::IO, mpr::MetaPropResult{:fixed})
    println(io, "  Meta-proportion fixed-effect result:")
    println(io, "  Weights (%): $(round.(mpr.wts ./ sum(mpr.wts) .* 100, sigdigits = 5))")
    println(io, "  Estimate: $(round(mpr.est, sigdigits = 6))")
    println(io, "  Variance (Std. error): $(round(mpr.var, sigdigits = 6)) ($(round(sqrt(mpr.var), sigdigits = 6)))")
    if mpr.data.metric in (:or, :rr)
    println(io, "  Exp(Estimate): $(round(exp(mpr.est), sigdigits = 6))")
    end
    println(io, "  Chi²: $(round(mpr.chisq, sigdigits = 6))")
    print(io,   "  Q: $(round(mpr.hetq, sigdigits = 6))")
    #print(io, "  I²: $(mpr.heti, sigdigits = 6))")
end

function Base.show(io::IO, mpr::MetaPropResult{:random})
    println(io, "  Meta-proportion random-effect result:")
    println(io, "  Weights (%): $(round.(mpr.wts ./ sum(mpr.wts) .* 100, sigdigits = 5))")
    println(io, "  Estimate: $(round(mpr.est, sigdigits = 6))")
    println(io, "  Variance  (Std. error): $(round(mpr.var, sigdigits = 6)) ($(round(sqrt(mpr.var), sigdigits = 6)))")
    if mpr.data.metric in (:or, :rr)
    println(io, "  Exp(Estimate): $(round(exp(mpr.est), sigdigits = 6))")
    end
    println(io, "  Chi²: $(round(mpr.chisq, sigdigits = 6))")
    println(io, "  Q: $(round(mpr.hetq, sigdigits = 6))")
    println(io, "  I²: $(round(mpr.heti, sigdigits = 6))")
    print(io,   "  τ²: $(round(mpr.hettau, sigdigits = 6))")
end
