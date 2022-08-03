

"""
    diffci(x1, n1, x2, n2; level = 0.95, method = :default)

Proportion difference (x1 / n1 - x2 / n2) confidence interval.

* 'method'

- `:mn` | `:default`- Miettinen & Nurminen;  Miettinen, O. and Nurminen, M. (1985), Comparative analysis of two rates. Statist. Med., 4: 213-226. doi:10.1002/sim.4780040211;
- `:fm` | `:mee` - Mee maximum likelihood method; Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176
- `:wald` - Wald CI without CC;
- `:waldcc` - Wald CI with CC;
- `:nhs`  - Newcombes Hybrid (wilson) Score interval; Newcombe RG (1998), Interval Estimation for the Difference Between Independent Proportions: Comparison of Eleven Methods. Statistics in Medicine 17, 873-890;
- `:nhscc` - Newcombes Hybrid Score CC; Newcombe (1998);
- `:ac` -  Agresti-Caffo interval; Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000);
- `:ha` - Hauck-Andersen; Hauck, W. W., & Anderson, S. (1986). A Comparison of Large-Sample Confidence Interval Methods for the Difference of Two Binomial Probabilities. The American Statistician, 40(4), 318–322. doi:10.1080/00031305.1986.10475426 ;
- `:mover` - Method of variance estimates recovery;
- `:jeffrey` - Brown, Li's Jeffreys.

# Reference:

* Brown, L.D., Cai, T.T., and DasGupta, A. Interval estimation for a binomial proportion. Statistical Science, 16(2):101–117, 2001.
* Farrington, C. P. and Manning, G. (1990), “Test Statistics and Sample Size Formulae for Comparative Binomial Trials with Null Hypothesis of Non-zero Risk Difference or Non-unity Relative Risk,” Statistics in Medicine, 9, 1447–1454
* Li HQ, Tang ML, Wong WK. Confidence intervals for ratio of two Poisson rates using the methodof variance estimates recovery. Computational Statistics 2014; 29(3-4):869-889
* Brown, L., Cai, T., & DasGupta, A. (2003). INTERVAL ESTIMATION IN EXPONENTIAL FAMILIES. Statistica Sinica, 13(1), 19-49.
"""
function diffci(x1, n1, x2, n2; level = 0.95, method = :default)

    alpha    = 1 - level
    if method == :mn || method == :default
        ci_diff_mn(x1, n1, x2, n2, alpha)
    elseif method == :fm || method == :mee
        ci_diff_fm(x1, n1, x2, n2, alpha)
    elseif method == :wald
        ci_diff_wald(x1, n1, x2, n2, alpha)
    elseif method == :waldcc
        ci_diff_wald_cc(x1, n1, x2, n2, alpha)
    elseif method == :nhs
        ci_diff_nhs(x1, n1, x2, n2, alpha)
    elseif method == :nhscc
        ci_diff_nhs_cc(x1, n1, x2, n2, alpha)
    elseif method == :ac
        ci_diff_ac(x1, n1, x2, n2, alpha)
    elseif method == :ha
        ci_diff_ha(x1, n1, x2, n2, alpha)
    elseif method == :mover
        ci_diff_mover(x1, n1, x2, n2, alpha)
    elseif method == :jeffrey
        ci_diff_jeffrey(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
"""
    diffci(contab::ConTab; level = 0.95, method = :default)

ConTab 2X2:

```
A | B
--|--
C | D
```

Difference: A / (A + B) - C / (C + D)
"""
function diffci(contab::ConTab; level = 0.95, method = :default)
    if !(size(contab.tab, 1) == size(contab.tab, 2) == 2) throw(ArgumentError("CI only for 2 X 2 tables.")) end
    alpha    = 1 - level
    x1 = contab.tab[1,1]
    n1 = x1 + contab.tab[1,2]
    x2 = contab.tab[2,1]
    n2 = x2 + contab.tab[2,2]
    diffci(x1, n1, x2, n2; level = level, method = method)
end
"""
    orci(x1, n1, x2, n2; level = 0.95, method = :default)

- `:mn` - MN Score (Miettinen&Nurminen, 1985);
- `:fm` | `:mee` - FM (same as MN Score, but not multiplied on `(n1 + n2) * (n1 + n2 - 1)`) (Mee RW, 1984; Farrington&Manning, 1990);
- `:woolf` - Woolf logit (Woolf, 1955);
- `:awoolf` - Adjusted Woolf interval (Gart adjusted logit) (Gart, 1966; Lawson, 2005);
- `:mover` - Method of variance estimates recovery (MOVER) (Donner&Zou, 2012);

Reference:

* Miettinen O. S., Nurminen M. (1985) Comparative analysis of two rates.Statistics in Medicine 4,213–226;
* Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics 40:1175-1176;
* Farrington, C. P. and Manning, G. (1990), “Test Statistics and Sample Size Formulae for Comparative Binomial Trials with Null Hypothesis of Non-zero Risk Difference or Non-unity Relative Risk,” Statistics in Medicine, 9, 1447–1454;
* Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253;
* Gart, J. J. (1966). Alternative analyses of contingency tables. Journal of the Royal Statistical Society. Series B (Methodological), 28:164-179;
* Lawson, R (2005). Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation and Computation, 33, 1095-1113;
* Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.

"""
function orci(x1, n1, x2, n2; level = 0.95, method = :default)
    alpha    = 1 - level
    if method == :mn || method == :default
        ci_or_mn(x1, n1, x2, n2, alpha)
    elseif method == :fm || method == :mee
        ci_or_fm(x1, n1, x2, n2, alpha)
    elseif method == :woolf
        ci_or_woolf(x1, n1, x2, n2, alpha)
    elseif method == :awoolf
        ci_or_awoolf(x1, n1, x2, n2, alpha)
    elseif method == :mover
        ci_or_mover(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
"""
    orci(contab::ConTab; level = 0.95, method = :default)

```
A | B
--|--
C | D
```

Odd ratio: (A / B) / (C / D)

"""
function orci(contab::ConTab; level = 0.95, method = :default)
    if !(size(contab.tab, 1) == size(contab.tab, 2) == 2) throw(ArgumentError("CI only for 2 X 2 tables.")) end
    x1 = contab.tab[1, 1]
    n1 = x1 + contab.tab[1, 2]
    x2 = contab.tab[2, 1]
    n2 = x2 + contab.tab[2, 2]
    orci(x1, n1, x2, n2; level = level, method = method)
end
"""
    rrci(x1, n1, x2, n2; level = 0.95, method = :default)

- `:mn` - Miettinen-Nurminen Score interval (Miettinen&Nurminen, 1985);
- `:fm` | `:mee` - FM Score interval (Mee RW, 1984; Farrington&Manning, 1990);
- `:cli` - Crude log interval, Gart (Gart&Nam, 1988);
- `:li` | `:wald` - Log interval / Katz / Wald interval (Katz et al, 1978);
- `:mover` - Method of variance estimates recovery (Donner&Zou, 2012);

Reference:

* Miettinen, O. and Nurminen, M. (1985), Comparative analysis of two rates. Statist. Med., 4: 213-226. doi:10.1002/sim.4780040211;
* Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics 40:1175-1176;
* Farrington, C. P., & Manning, G. (1990). Test statistics and sample size formulae for comparative binomial trials with null hypothesis of non-zero risk difference or non-unity relative risk. Statistics in Medicine, 9(12), 1447–1454. doi:10.1002/sim.4780091208;
* Gart, JJ and Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338;
* Katz D, Baptista J, Azen SP and Pike MC. Obtaining confidence intervals for the risk ratio in cohort studies. Biometrics 1978; 34: 469–474;
* Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.
"""
function rrci(x1, n1, x2, n2; level = 0.95, method = :default)
    alpha    = 1 - level
    if method == :mn || method == :default
        ci_rr_mn(x1, n1, x2, n2, alpha)
    elseif method == :fm || method == :mee
        ci_rr_fm(x1, n1, x2, n2, alpha)
    elseif method == :cli
        ci_rr_cli(x1, n1, x2, n2, alpha)
    elseif method == :li || method == :wald
        ci_rr_li(x1, n1, x2, n2, alpha)
    elseif method == :mover
        ci_rr_mover(x1, n1, x2, n2, alpha)
    else
        throw(ArgumentError("unknown ci method=$(method)"))
    end
end
"""
    rrci(contab::ConTab; level = 0.95, method = :default)

```
A | B
--|--
C | D
```

Risk ratio: (A / (A + B)) / (C / (C + D)
"""
function rrci(contab::ConTab; level = 0.95, method = :default)
    if !(size(contab.tab, 1) == size(contab.tab, 2) == 2) throw(ArgumentError("CI only for 2 X 2 tables.")) end
    alpha    = 1 - level
    x1 = contab.tab[1, 1]
    n1 = x1 + contab.tab[1, 2]
    x2 = contab.tab[2, 1]
    n2 = x2 + contab.tab[2, 2]
    rrci(x1, n1, x2, n2; level = level, method = method)
end

"""
    propci(x::Int, n::Int; level = 0.95, method = :default)

`method`:

- `:wilson` | `:default` - Wilson's confidence interval (CI) for a single proportion (wilson score) (Wilson, 1927);
- `:wilsoncc` - Wilson's CI with continuity correction (CC);
- `:cp` - Clopper-Pearson exact CI (Clopper&Pearson, 1934);
- `:blaker` - Blaker exact CI for discrete distributions (Blaker, 2000);
- `:soc` - SOC: Second-Order corrected CI;
- `:arc` - Arcsine CI;
- `:wald` - Wald CI without CC;
- `:waldcc` - Wald CI with CC;
- `:ac` - Agresti-Coull;
- `:jeffrey` - Jeffreys interval.

Reference:

* Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inference J. Amer.Stat. Assoc 22, 209–212;
* Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the caseof the binomial.Biometrika26, 404–413;
* Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions, Canadian Journal of Statistics 28 (4), 783–798;

"""
function propci(x::Int, n::Int; level = 0.95, method = :default)
    if  x > n throw(ArgumentError("x > n")) end
    alpha    = 1 - level
    if method == :wilson || method == :default
        fx = ci_prop_wilson
    elseif method==:wilsoncc
        fx = ci_prop_wilson_cc
    elseif method==:cp
        fx = ci_prop_cp
    elseif method==:blaker
        fx = ci_prop_blaker
    elseif method==:soc
        fx = ci_prop_soc
    elseif method==:arc
        fx = ci_prop_arc
    elseif method==:wald
        fx = ci_prop_wald
    elseif method==:waldcc
        fx = ci_prop_wald_cc
    elseif method==:ac
        fx = ci_prop_ac
    elseif method==:jeffrey
        fx = ci_prop_jeffrey
    else
        throw(ArgumentError("unknown method!"))
    end
    fx(x, n, alpha)
end
"""
    propci(contab::ConTab; level = 0.95, method = :default)
"""
function propci(contab::ConTab; level = 0.95, method = :default)
    alpha    = 1 - level
    if  size(contab.tab, 2) != 2 throw(ArgumentError("CI only for N X 2 tables.")) end
    if size(contab, 1) > 1
        v = Vector{Tuple{Float64, Float64}}(undef, size(contab, 1))
        for i = 1:size(contab, 1)
            x = contab.tab[i, 1]
            n = x + contab.tab[i, 2]
            #println(x, " : ", n)
            v[i] = propci(x, n; level = level, method = method)
        end
        return v
    else
        x = contab.tab[1, 1]
        n = x + contab.tab[1, 2]
        return propci(x, n; level = level, method = method)
    end
end

"""
    propci(prop::Proportion; level = 0.95, method = :default)
"""
function propci(prop::Proportion; level = 0.95, method = :default)
    propci(prop.x, prop.n; level = level, method = method)
end

"""
    mpropci(contab::ConTab; level = 0.95, method = :default)

Multinomial proportions confidence interval.

`method`:

- `goodman` | `default` Goodman, L.A. (1965). On Simultaneous Confidence Intervals for Multinomial Proportions. Technometrics 7: 247-254.

"""
function mpropci(contab::ConTab; level = 0.95, method = :default)
    alpha    = 1 - level
    if  size(contab.tab, 2) <= 2 throw(ArgumentError("CI only for N X M tables where M > 2")) end

    if method == :goodman || method == :default
        fx = ci_prop_goodman
    else
        throw(ArgumentError("unknown method!"))
    end
    if size(contab, 1) > 1
        v = Vector{Vector{Tuple{Float64, Float64}}}(undef, size(contab, 1))
        for i = 1:size(contab, 1)
            v[i] = fx(view(contab.tab, i, :), alpha)
        end
        return v
    else
        return fx(view(contab.tab, 1, :), alpha)
    end
end

################################################################################
# Proportion difference CI
################################################################################
# Method of Mee 1984 with Miettinen and Nurminen modification n / (n - 1) Newcombe 1998
# Score intervals for the difference of two binomial proportions
@inline function mle_diff(p1, n1, p2, n2, δ)
    if p1 - p2 - δ == 0 return 0.0 end
    θ = n2 / n1
    a = 1 + θ
    b = -(1 + θ + p1 + θ * p2 + δ * (θ + 2))
    c = δ^2 + δ * (2p1 + θ + 1) + p1 + θ * p2
    d = -p1 * δ * (1 + δ)
    v = (b / 3a)^3 - b * c / (6 * a * a) + d / 2a
    u = sign(v) * sqrt((b / 3a)^2 - c / 3a)
    w = (pi + acos(v / u^3)) / 3
    p1n = 2u * cos(w) - b / 3a
    p2n = p1n - δ
    return p1n, p2n
end
@inline function mn_fm_diff_z_val(p1, n1, p2, n2, est, δ)
    p1n, p2n = mle_diff(p1, n1, p2, n2, δ)
    return (est - δ)^2 / (p1n * (1 - p1n) / n1 + p2n * (1 - p2n) / n2)
end
@inline function mn_diff_z_val(p1, n1, p2, n2, est, δ)
    mn_fm_diff_z_val(p1, n1, p2, n2, est, δ) / (n1 + n2) * (n1 + n2 - 1)
end
function ci_diff_mn(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    lcis, ucis = ci_diff_nhs_cc(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z²        = quantile(Chisq(1), 1 - alpha)
    fmnd(x)  = mn_diff_z_val(p1, n1, p2, n2, est, x) - z²
    #=
    if fmnd(lcis) * fmnd(est - eps()) < 0.0
        ll = lcis
        lu = est - eps()
    else
        ll = -1.0 + eps()
        lu = lcis
    end
    if fmnd(ucis) * fmnd(est + eps()) < 0.0
        ul = est + eps()
        uu = ucis
    else
        ul = ucis
        uu = 1.0 - eps()
    end
    =#
    lci = find_zero(fmnd, lcis)
    uci = find_zero(fmnd, ucis)
    return lci, uci
end
# FM / MEE
# Mee RW (1984) Confidence bounds for the difference between two probabilities, Biometrics40:1175-1176
# MN - no correction
function ci_diff_fm(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    lcis, ucis = ci_diff_nhs_cc(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z²        = quantile(Chisq(1), 1 - alpha)
    fmnd(x)  = mn_fm_diff_z_val(p1, n1, p2, n2, est, x) - z²
    #=
    if fmnd(lcis) * fmnd(est - eps()) < 0.0
        ll = lcis
        lu = est - eps()
    else
        ll = -1.0 + eps()
        lu = lcis
    end
    if fmnd(ucis) * fmnd(est + eps()) < 0.0
        ul = est + eps()
        uu = ucis
    else
        ul = ucis
        uu = 1.0 - eps()
    end
    =#
    lci = find_zero(fmnd, lcis)
    uci = find_zero(fmnd, ucis)
    return  lci, uci
end
# Wald
function ci_diff_wald(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return est - z * se, est + z * se
end
# Wald continuity correction
function ci_diff_wald_cc(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    cc       = (1 / n1 + 1 / n2) / 2
    return est - z * se - cc, est + z * se + cc
end
# Newcombes Hybrid (wilson) Score interval for the difference of proportions
# Newcombe 1998
function ci_diff_nhs(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    lci1, uci1 = ci_prop_wilson(x1, n1, alpha)
    lci2, uci2 = ci_prop_wilson(x2, n2, alpha)
    return  est - z * sqrt(lci1 * (1 - lci1)/n1 + uci2 * (1 - uci2) / n2), est + z * sqrt(uci1 * (1 - uci1) / n1 + lci2 * (1 - lci2) / n2)
end
# Newcombes Hybrid Score continuity correction
function ci_diff_nhs_cc(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    lci1, uci1 = ci_prop_wilson_cc(x1, n1, alpha)
    lci2, uci2 = ci_prop_wilson_cc(x2, n2, alpha)
    return  est - sqrt((p1 - lci1)^2 + (uci2 - p2)^2), est + sqrt((uci1 - p1)^2 + (p2 - lci2)^2)
end
# Agresti-Caffo interval for the difference of proportions
# Agresti A, Caffo B., “Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures”, American Statistician 54: 280–288 (2000)
function ci_diff_ac(x1, n1, x2, n2, alpha)
    z        = quantile(Normal(), 1 - alpha / 2)
    p1I      = (x1 + 1) / (n1 + 2)
    p2I      = (x2 + 1) / (n2 + 2)
    n1I      = n1 + 2
    n2I      = n2 + 2
    est      = p1I - p2I
    se       = sqrt(p1I * (1 - p1I) / n1I + p2I * (1 - p2I) / n2I)
    return  est - z * se, est + z * se
end
# Hauck-Andersen
# Hauck, W. W., & Anderson, S. (1986). A Comparison of Large-Sample Confidence Interval Methods for the Difference of Two Binomial Probabilities. The American Statistician, 40(4), 318–322. doi:10.1080/00031305.1986.10475426

function ci_diff_ha(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    z        = quantile(Normal(), 1 - alpha / 2)
    se       = sqrt(p1 * (1 - p1) / (n1 - 1) + p2 * (1 - p2) / (n2 - 1))
    cc       = 1 / min(n1, n2)
    return est - z * se - cc, est + z * se + cc
end
# Brown, Li's Jeffreys
function ci_diff_jeffrey(x1, n1, x2, n2, alpha)
    p1   = (x1 +  1/2) / (n1 + 1)
    p2   = (x2 +  1/2) / (n2 + 1)
    se   = sqrt(p1*(1 - p1) / n1 + p2 * (1 - p2) / n2)
    z    = quantile(Normal(), 1 - alpha / 2)
    est  = p1 - p2
    return  max(-1.0, est - z * se), min(1.0, est + z * se)
end
# Method of variance estimates recovery
function ci_diff_mover(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 - p2
    Z        = quantile(Normal(), 1 - alpha / 2)
    lci1, uci1   = ci_prop_wilson(x1, n1, alpha)
    lci2, uci2   = ci_prop_wilson(x2, n2, alpha)
    lci      = est - sqrt((p1 - lci1)^2 + (uci2 - p2)^2)
    uci      = est + sqrt((uci1 - p1)^2 + (p2 - lci2)^2)
    return lci, uci
end
################################################################################
# Odd ratio CI
################################################################################
@inline function mle_or(φ, x1, n1, x2, n2)
    a  = n2 * (φ - 1)
    b  = φ * n1 + n2 - (x1 + x2) * (φ - 1)
    c  = -(x1 + x2)
    p2 = (-b + sqrt(b * b - 4 * a * c)) / 2a
    p1 = p2 * φ / (1 + p2 * (φ - 1))
    return p1, p2
end
@inline function mle_fm_or_z_val(φ, x1, n1, x2, n2)
    p1 = x1 / n1
    pmle1, pmle2 = mle_or(φ, x1, n1, x2, n2)
    return (n1 * (p1 - pmle1))^2 * (1 / (n1 * pmle1 * (1 - pmle1)) + 1/(n2 * pmle2 * (1 - pmle2)))
end
@inline function mle_or_z_val(φ, x1, n1, x2, n2)
    mle_fm_or_z_val(φ, x1, n1, x2, n2) / (n1 + n2) * (n1 + n2 - 1)
end
################################################################################
# Miettinen O. S., Nurminen M. (1985) Comparative analysis of two rates.Statistics in Medicine 4,213–226
# MN Score
function ci_or_mn(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    z²        = quantile(Chisq(1), 1 - alpha)
    fmnor(x) = mle_or_z_val(x, x1, n1, x2, n2) - z²
    lci, uci = ci_or_awoolf(x1, n1, x2, n2, alpha)
    if (x1== 0 && x2 == 0) || (x1 == n1 && x2 == n2)
        return 0.0, Inf
    elseif x1==0 || x2 == n2
        return 0.0, find_zero(fmnor, uci)
    elseif x1 == n1 || x2 == 0
        return find_zero(fmnor, lci), Inf
    else
        est  = (x1 / (n1 - x1)) / (x2 / (n2 - x2))
        return find_zero(fmnor, lci), find_zero(fmnor, uci)
    end
end
# FM Score
# Mee RW (1984) Confidence bounds for the difference between two probabilities,Biometrics40:1175-1176
# Farrington, C. P. and Manning, G. (1990), “Test Statistics and Sample Size Formulae for Comparative Binomial Trials with Null Hypothesis of Non-zero Risk Difference or Non-unity Relative Risk,” Statistics in Medicine, 9, 1447–1454

function ci_or_fm(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    z²        = quantile(Chisq(1), 1 - alpha)
    fmnor(x) = mle_fm_or_z_val(x, x1, n1, x2, n2) - z²
    lci, uci = ci_or_awoolf(x1, n1, x2, n2, alpha)
    if (x1== 0 && x2 == 0) || (x1 == n1 && x2 == n2)
        return 0.0, Inf
    elseif x1==0 || x2 == n2
        return 0.0, find_zero(fmnor, uci)
    elseif x1 == n1 || x2 == 0
        return find_zero(fmnor, lci), Inf
    else
        est  = (x1 / (n1 - x1)) / (x2 / (n2 - x2))
        return find_zero(fmnor, lci), find_zero(fmnor, uci)
    end
end
# Woolf logit
# Woolf, B. (1955). On estimating the relation between blood group and disease. Annals of human genetics, 19(4):251-253.
function ci_or_woolf(x1, n1, x2, n2, alpha)
    xa        = x1
    xb        = n1 - x1
    xc        = x2
    xd        = n2 - x2
    est       = xa*xd/xc/xb
    estI      = log(est)
    se        = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
    z         = quantile(Normal(), 1 - alpha / 2)
    return exp(estI - z * se), exp(estI + z * se)
end
# Adjusted Woolf interval (Gart adjusted logit) Lawson, R (2005):Smallsample confidence intervals for the odds ratio.  Communication in Statistics Simulation andComputation, 33, 1095-1113.
# Gart, J. J. (1966). Alternative analyses of contingency tables. Journal of the Royal Statistical Society. Series B (Methodological), 28:164-179.
function ci_or_awoolf(x1, n1, x2, n2, alpha)
    xa        = x1 +  1/2
    xb        = n1 - x1 +  1/2
    xc        = x2 + 1/2
    xd        = n2 - x2 +  1/2
    est       = xa*xd/xc/xb
    estI      = log(est)
    se        = sqrt(1/xa + 1/xb + 1/xc + 1/xd)
    z         = quantile(Normal(), 1 - alpha/2)
    return exp(estI - z * se), exp(estI + z * se)
end
# Method of variance estimates recovery
# Donner, A. and Zou, G. (2012). Closed-form confidence intervals for functions of the normal mean and standard deviation. Statistical Methods in Medical Research, 21(4):347-359.
function ci_or_mover(x1, n1, x2, n2, alpha)
    p1       = (x1/(n1-x1))
    p2       = (x2/(n2-x2))
    est      = p1/p2
    z        = quantile(Normal(), 1 - alpha / 2)
    lci1, uci1   = ci_prop_wilson(x1, n1, alpha)
    lci2, uci2   = ci_prop_wilson(x2, n2, alpha)
    vl1      = lci1/(1 - lci1)
    vu1      = uci1/(1 - uci1)
    vl2      = lci2/(1 - lci2)
    vu2      = uci2/(1 - uci2)
    lci      = (p1*p2-sqrt((p1*p2)^2 - vl1*vu2*(2*p1-vl1)*(2*p2-vu2)))/(vu2*(2*p2 - vu2))
    uci      = (p1*p2+sqrt((p1*p2)^2 - vu1*vl2*(2*p1-vu1)*(2*p2-vl2)))/(vl2*(2*p2 - vl2))
    return lci, uci
end
################################################################################
# Risk ratio CI
################################################################################
@inline function mle_rr(φ, x1, n1, x2, n2)
    a = (n1 + n2) * φ
    b = -(φ * (x2 + n1) + x1 + n2)
    c = x1 + x2
    p2 = (-b - sqrt(b * b - 4a * c)) / 2a
    p1 = p2 * φ
    return p1, p2
end
@inline function mle_fm_rr_z_val(φ, x1, n1, x2, n2)
    p1 = x1 / n1
    p2 = x2 / n2
    pmle1, pmle2 = mle_rr(φ, x1, n1, x2, n2)
    return (p1 - φ * p2)^2 / (pmle1 * (1 - pmle1) / n1 + φ * φ * pmle2 * (1 - pmle2) / n2)
end
@inline function mle_rr_z_val(φ, x1, n1, x2, n2)
    mle_fm_rr_z_val(φ, x1, n1, x2, n2) / (n1 + n2) * (n1 + n2 - 1)
end
################################################################################
# Miettinen-Nurminen Score interval
# Miettinen, O. and Nurminen, M. (1985), Comparative analysis of two rates. Statist. Med., 4: 213-226. doi:10.1002/sim.4780040211
function ci_rr_mn(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    lci, uci = ci_rr_cli(x1, n1, x2, n2, alpha, 1/4)
    z²        = quantile(Chisq(1), 1 - alpha)
    fmnrr(x) = mle_rr_z_val(x, x1, n1, x2, n2) - z²
    if (x1 == 0 && x2 == 0) || (x1==n1 && x2==n2)
        return  0.0, Inf
    elseif x1 == 0 || x2 == n2
        return 0.0, find_zero(fmnrr, uci)
    elseif x1 == n1 || x2 == 0
        return find_zero(fmnrr, lci), Inf
    else
        est = (x1 / n1) / (x2 / n2)
        return find_zero(fmnrr, lci), find_zero(fmnrr, uci)
    end
end
# FM Score interval
# Farrington, C. P., & Manning, G. (1990). Test statistics and sample size formulae for comparative binomial trials with null hypothesis of non-zero risk difference or non-unity relative risk. Statistics in Medicine, 9(12), 1447–1454. doi:10.1002/sim.4780091208
function ci_rr_fm(x1, n1, x2, n2, alpha; atol::Float64 = 1E-8)
    lci, uci = ci_rr_cli(x1, n1, x2, n2, alpha, 1/4)
    z²           = quantile(Chisq(1), 1 - alpha)
    fmnrr(x)     = mle_fm_rr_z_val(x, x1, n1, x2, n2) - z²
    if (x1 == 0 && x2 == 0) || (x1 == n1 && x2 == n2)
        return  0.0, Inf
    elseif x1 == 0 || x2 == n2
        return 0.0, find_zero(fmnrr, uci)
    elseif x1 == n1 || x2 == 0
        return find_zero(fmnrr, lci), Inf
    else
        #est = (x1 / n1) / (x2 / n2)
        return find_zero(fmnrr, lci), find_zero(fmnrr, uci)
    end
end
# Crude log interval
# Gart, JJand Nam, J (1988): Approximate interval estimation of the ratio of binomial parameters: Areview and corrections for skewness. Biometrics 44, 323-338.
function ci_rr_cli(x1, n1, x2, n2, alpha, adj = 1/2)
    x1I       = x1 + adj
    x2I       = x2 + adj
    n1I       = n1 + adj
    n2I       = n2 + adj
    estI      = log((x1I / n1I) / (x2I / n2I))
    se        = sqrt(1 / x2I + 1 / x1I - 1 / n2I - 1 / n1I)
    est       = (x1 / n1) / (x2 / n2)
    z         =  quantile(Normal(), 1 - alpha / 2)
    return exp(estI - z * se), exp(estI + z * se)
end
# Katz D, Baptista J, Azen SP and Pike MC. Obtaining confidence intervals for the risk ratio in cohort studies. Biometrics 1978; 34: 469–474
function ci_rr_li(x1, n1, x2, n2, alpha)
    est       = (x1 / n1) / (x2 / n2)
    estI      = log(est)
    se        = sqrt(1 / x2 + 1 / x1 - 1 / n2 - 1 / n1)
    z         = quantile(Normal(), 1 - alpha / 2)
    return exp(estI - z * se), exp(estI + z * se)
end
# Method of variance estimates recovery (Donner, Zou, 2012)
function ci_rr_mover(x1, n1, x2, n2, alpha)
    p1       = x1 / n1
    p2       = x2 / n2
    est      = p1 / p2
    Z        = quantile(Normal(), 1 - alpha / 2)
    lci1, uci1   = ci_prop_wilson(x1, n1, alpha)
    lci2, uci2   = ci_prop_wilson(x2, n2, alpha)
    lci      = (p1 * p2 - sqrt((p1 * p2)^2 - lci1 * uci2 * (2 * p1 - lci1) * (2 * p2 - uci2))) / (uci2 * (2 * p2 - uci2))
    uci      = (p1 * p2 + sqrt((p1 * p2)^2 - uci1 * lci2 * (2 * p1 - uci1) * (2 * p2 - lci2))) / (lci2 * (2 * p2 - lci2))
    return  lci, uci
end
################################################################################
# Prorotion CI
################################################################################
# Wilson’s confidence interval for a single proportion, wilson score
# Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inferenceJ. Amer.Stat. Assoc22, 209–212
function ci_prop_wilson(x, n, alpha)
    z   = abs(quantile(Normal(), 1 - alpha / 2))
    p   = x / n
    d   = 1 + (z^2) / n
    se  = sqrt((p * (1 - p) + (z^2) / (4 * n)) / n) / d
    est = (p + (z^2) / (2 * n)) / d
    return est - z * se, est + z * se
end
# Wilson CC
# Newcombe, R. G. (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods". Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 959561
function ci_prop_wilson_cc(x, n, alpha)
    z = abs(quantile(Normal(), 1 - alpha / 2))
    p = x / n
    l = (2*n*p+z*z-1-z*sqrt(z*z-2-1/n+4*p*(n*(1-p)+1)))/2/(n+z*z)
    u = (2*n*p+z*z+1+z*sqrt(z*z+2-1/n+4*p*(n*(1-p)-1)))/2/(n+z*z)
    return min(p, l), max(p, u)
end
#Clopper-Pearson exatct CI
#Clopper, C. and Pearson, E.S. (1934) The use of confidence or fiducial limits illustrated in the caseof the binomial.Biometrika26, 404–413.
function ci_prop_cp(x, n, alpha)
    if x == 0
        ll = 0.0
        ul = 1 - (alpha / 2)^(1 / n)
    elseif x == n
        ul = 1.0
        ll = (alpha / 2)^(1 / n)
    else
        ll = 1/(1 + (n - x + 1)/(x * quantile(FDist(2 * x, 2 * (n - x + 1)), alpha / 2)))
        ul = 1/(1 + (n - x) / ((x + 1) * quantile(FDist(2 * (x + 1), 2 * (n - x)), 1 - alpha / 2)))
    end
    return ll, ul
end
# Blaker CI
# Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions,Canadian Journal of Statistics28 (4), 783–798
function ci_prop_blaker(x, n, alpha; atol::Float64 = 1E-8)
    lower = 0.; upper = 1.;
    fx(p) =  acceptbin(x, n, p) - alpha
    if n != 0
        lower = quantile(Beta(x, n - x + 1), alpha / 2)
        lower = find_zero(fx, lower)
    end
    if x != n
        upper = quantile(Beta(x + 1, n - x), 1 - alpha / 2)
        upper = find_zero(fx, upper)
    end
    return lower,upper
end
@inline function acceptbin(x, n, p)
    BIN = Binomial(n,p)
    p1 = 1-cdf(BIN,x-1)
    p2 =   cdf(BIN,x)
    a1 = p1 + cdf(BIN, quantile(BIN,p1)-1)
    a2 = p2+1-cdf(BIN, quantile(BIN,1-p2))
    return min(a1,a2)
end
# SOC  Second-Order corrected
# T. Tony Cai One-sided confdence intervals in discrete distributions doi:10.1016/j.jspi.2004.01.00
function ci_prop_soc(x, n, alpha)
    p  = x / n
    k  = quantile(Normal(), 1 - alpha / 2)
    k2 = k^2
    η  = k2 / 3 + 1 / 6
    γ1 = -(k2 * 13 / 18 + 17 / 18)
    γ2 = k2 / 18 + 7 / 36
    m  = (x + η) / (n + 2 * η)
    b  = k * sqrt(p * (1 - p) + (γ1 * p *(1 - p) + γ2) / n) / sqrt(n)
    return m - b, m + b
end
# Arcsine
function ci_prop_arc(x, n, alpha)
    q = quantile(Normal(), 1-alpha/2)
    p = x/n
    z = q/(2*sqrt(n))
    return sin(asin(sqrt(p))-z)^2, sin(asin(sqrt(p))+z)^2
end
# Wald CI
function ci_prop_wald(x, n, alpha)
    p=x/n
    b = quantile(Normal(), 1-alpha/2)*sqrt(p*(1-p)/n)
    return p-b, p+b
end
# Wald CI CC
function ci_prop_wald_cc(x::Int, n::Int, alpha::Real)
    p=x/n
    b = quantile(Normal(), 1-alpha/2)*sqrt(p*(1-p)/n)
    cc = 0.5 / n
    return p-b-cc, p+b+cc
end
# Agresti-Coull
function ci_prop_ac(x, n, alpha)
    z = quantile(Normal(), 1 - alpha / 2)
    n = n + z^2
    est = (x + z ^ 2 / 2) / n
    se = sqrt(est * (1 - est) / n)
    return (est - z * se, est + z * se)
end

# Jeffreys interval
function ci_prop_jeffrey(x, n, alpha)
    (quantile(Beta(x + 1/2, n - x + 1/2), alpha/2), quantile(Beta(x + 1/2, n - x + 1/2), 1-alpha/2))
end

################################################################################
# Goodman CI
################################################################################
# Goodman, L.A. (1965). On Simultaneous Confidence Intervals for Multinomial Proportions. Technometrics 7: 247-254.
# https://blogs.sas.com/content/iml/2017/02/15/confidence-intervals-multinomial-proportions.html
# https://rdrr.io/cran/CoinMinD/man/GM.html
function ci_prop_goodman(v, alpha::T2) where T2
    k   = length(v)
    s   = sum(v)
    p   = v ./ s
    chi = quantile(Chisq(one(Int)), one(T2) - alpha/k)
    ci  = Vector{Tuple{Float64, Float64}}(undef, k)
    d   = 2 * (chi + s)
     @inbounds @simd for i = 1:k
        ci[i] = ( (chi + 2v[i] - sqrt(chi * chi + 4v[i] * chi*(1.0 - v[i] / s))) / d , (chi + 2v[i] + sqrt(chi * chi + 4v[i] * chi * (1.0 - v[i] / s))) / d )
    end
    ci
end
