
using Test
using DataFrames, CSV, CategoricalArrays, HypothesisTests

path     = dirname(@__FILE__)
io       = IOBuffer();


function test_ci_approx(x, y)
        all(map(isapprox, x, y))
end


# library(DescTools) for R using to get reference values
#=
int = "clopper-pearson"
BinomCI(26,78, method = int)
BinomCI(38,100, method = int)
BinomCI(0,5, method = int)
=#

#dat.li2007
#dat.hine1989
#dat.graves2010
#dat.bourassa1996
#https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf

@testset "  Proportion Confidence Intarvals                          " begin
    ct = MetidaFreq.contab([38, 62])
    t   = BinomialTest(26, 78)
    t2  = BinomialTest(0, 5)

    # Wald
#=
> BinomCI(26,78, method = "wald")
           est    lwr.ci    upr.ci
[1,] 0.3333333 0.2287182 0.4379485
> BinomCI(38,100, method = "wald")
      est   lwr.ci   upr.ci
[1,] 0.38 0.284866 0.475134
> BinomCI(0,5, method = "wald")
     est lwr.ci upr.ci
[1,]   0      0      0
=#
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [0.28486600512143206, 0.47513399487856794] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [0.28486600512143206, 0.47513399487856794] atol=1E-6
    @test test_ci_approx(confint(t,  method=:wald), MetidaFreq.propci(26, 78; method = :wald))
    @test test_ci_approx(confint(t2, method=:wald), MetidaFreq.propci(0, 5; method = :wald))

    # Wald cc
    #=
> BinomCI(26,78, method = "waldcc")
           est    lwr.ci    upr.ci
[1,] 0.3333333 0.2223079 0.4443587
> BinomCI(38,100, method = "waldcc")
      est   lwr.ci   upr.ci
[1,] 0.38 0.279866 0.480134
> BinomCI(0,5, method = "waldcc")
     est lwr.ci upr.ci
[1,]   0      0    0.1
    =#
    ci = MetidaFreq.propci(ct; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [0.27986600512143206, 0.48013399487856795] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [0.27986600512143206, 0.48013399487856795] atol=1E-6
    ci = MetidaFreq.propci(0, 5; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [0, 0.1] atol=1E-6 # Shold be fixed?

    # Wilson
#=
> int = "wilson"
> BinomCI(26,78, method = int)
           est    lwr.ci    upr.ci
[1,] 0.3333333 0.2387267 0.4435859
> BinomCI(38,100, method = int)
      est   lwr.ci    upr.ci
[1,] 0.38 0.290976 0.4779024
> BinomCI(0,5, method = int)
     est lwr.ci    upr.ci
[1,]   0      0 0.4344825
> 
=#    
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wilson)
    @test collect(ci)  ≈ [0.2909759925247873, 0.47790244704488943] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wilson)
    @test collect(ci)  ≈ [0.2909759925247873, 0.47790244704488943] atol=1E-6
    @test test_ci_approx(confint(t,  method=:wilson), MetidaFreq.propci(26, 78; method = :wilson))
    @test test_ci_approx(confint(t2, method=:wilson), MetidaFreq.propci(0, 5; method = :wilson))

#=
> int = "wilsoncc"
> BinomCI(26,78, method = int)
           est   lwr.ci    upr.ci
[1,] 0.3333333 0.233094 0.4501519
> BinomCI(38,100, method = int)
      est    lwr.ci    upr.ci
[1,] 0.38 0.2863947 0.4829411
> BinomCI(0,5, method = int)
     est lwr.ci   upr.ci
[1,]   0      0 0.537056
=#
    # Wilson cc - check
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wilsoncc)
    @test collect(ci)  ≈ [0.28639471634406627, 0.48294114679105093] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wilsoncc)
    @test collect(ci)  ≈ [0.28639471634406627, 0.48294114679105093] atol=1E-6
    ci = MetidaFreq.propci(26, 78; method = :wilsoncc)
    @test collect(ci)  ≈ [0.2330940074016006, 0.4501518829756797] atol=1E-6
    ci = MetidaFreq.propci(0, 5; method = :wilsoncc)
    @test collect(ci)  ≈ [0.0, 0.537056017467524] atol=1E-6
    ci = MetidaFreq.propci(5, 5; method = :wilsoncc)
    @test collect(ci)  ≈ [0.46294398253247615, 1] atol=1E-6

 #=
 > int = "clopper-pearson"
> BinomCI(26,78, method = int)
           est    lwr.ci    upr.ci
[1,] 0.3333333 0.2305852 0.4491667
> BinomCI(38,100, method = int)
      est    lwr.ci    upr.ci
[1,] 0.38 0.2847675 0.4825393
> BinomCI(0,5, method = int)
     est lwr.ci    upr.ci
[1,]   0      0 0.5218238
> 
 =#   
    # Clopper-Pirson
    ci = MetidaFreq.propci(ct; level = 0.95, method = :cp)
    @test collect(ci)  ≈ [0.2847674761414794, 0.48253930575080645] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :cp)
    @test collect(ci)  ≈ [0.2847674761414794, 0.48253930575080645] atol=1E-6

    @test_nowarn MetidaFreq.propci(0, 100; level = 0.95, method = :cp)
    @test_nowarn MetidaFreq.propci(100, 100; level = 0.95, method = :cp)


    # SOC
    ci = MetidaFreq.propci(ct; level = 0.95, method = :soc)
    @test collect(ci)  ≈ [0.2891917018839231, 0.4775592393403469] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :soc)
    @test collect(ci)  ≈ [0.2891917018839231, 0.4775592393403469] atol=1E-6

    # Blaker
    # Validate
    ci = MetidaFreq.propci(ct; level = 0.95, method = :blaker)
    @test collect(ci)  ≈ [0.2881924139040972, 0.4798218835425426] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :blaker)
    @test collect(ci)  ≈ [0.2881924139040972, 0.4798218835425426] atol=1E-6

    # arcsine
    ci = MetidaFreq.propci(ct; level = 0.95, method = :arc)
    @test collect(ci)  ≈ [0.2877714314998773, 0.47682358116201534] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :arc)
    @test collect(ci)  ≈ [0.2877714314998773, 0.47682358116201534] atol=1E-6

    # Jeffrey
    ci = MetidaFreq.propci(ct; level = 0.95, method = :jeffrey)
    @test collect(ci)  ≈ [0.2893837310094326, 0.4774506923359312] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :jeffrey)
    @test collect(ci)  ≈ [0.2893837310094326, 0.4774506923359312] atol=1E-6

    # AC
    # Validate
    ci = MetidaFreq.propci(ct; level = 0.95, method = :ac)
    @test collect(ci)  ≈ [0.2908745228985889, 0.4780039166710877] atol=1E-6

    ############################################################################
    # default = wilson, level = 0.95
    ci = MetidaFreq.propci(MetidaFreq.Proportion(38, 100))
    @test collect(ci)  ≈ [0.2909759925247873, 0.47790244704488943] atol=1E-6
end

@testset "  Proportion difference Confidence Intarvals               " begin
    ct = MetidaFreq.contab([30 70; 40 50])

    # Miettinen & Nurminen
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :mn)
    @test collect(ci)  ≈ [-0.2781290897168457, -0.006708341755865329] atol=1E-6
    @test_nowarn MetidaFreq.diffci(0, 100, 1, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.diffci(1, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.diffci(0, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.diffci(100, 100, 100, 100; level = 0.95, method = :mn)

    # FM | Mee
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [-0.27778778468650384, -0.007071207814437397] atol=1E-6
    @test_nowarn MetidaFreq.diffci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.diffci(1, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.diffci(0, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.diffci(100, 100, 100, 100; level = 0.95, method = :fm)

    # Wald
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [-0.28084842238, -0.00804046650] atol=1E-6

    # Wald CC
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [-0.29140397794, 0.00251508905] atol=1E-6

    # Newcombe Hybrid Score
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :nhs)
    @test collect(ci)  ≈ [-0.275381800, -0.007158419] atol=1E-6

    # Newcombe Hybrid Score CC
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :nhscc)
    @test collect(ci)  ≈ [-0.2823838433576864, 0.00020433786401141685] atol=1E-6

    # Agresti-Caffo interval
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :ac)
    @test collect(ci)  ≈ [-0.276944506, -0.006516705] atol=1E-6

    # ha
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :ha)
    @test collect(ci)  ≈ [-0.2926903294586995, 0.003801440569810589] atol=1E-6

    # mover
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :mover)
    @test collect(ci)  ≈ [-0.2753818003977219, -0.007158418963689267] atol=1E-6

    # jeffrey
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :jeffrey)
    @test collect(ci)  ≈ [-0.27960020759495857, -0.006549286475327543] atol=1E-6


end

@testset "  Odd ratio Confidence Intarvals                           " begin
    ct = MetidaFreq.contab([30 70; 40 50])

    # Miettinen & Nurminen
    # 0.2953741424 - 0.9716669781
    ci = MetidaFreq.orci(ct; level = 0.95, method = :mn)
    @test collect(ci)  ≈ [0.2953741424, 0.9716669781] atol=1E-6
    @test_nowarn MetidaFreq.orci(0, 100, 1, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.orci(1, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.orci(0, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.orci(100, 100, 100, 100; level = 0.95, method = :mn)
    # FM
    # 0.2958336891 - 0.9701570484
    ci = MetidaFreq.orci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [0.2958336891, 0.9701570484] atol=1E-6
    @test_nowarn MetidaFreq.orci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.orci(1, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.orci(0, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.orci(100, 100, 100, 100; level = 0.95, method = :fm)
    # woolf | Wald
    # 0.2950420027 - 0.9727082695
    ci = MetidaFreq.orci(ct; level = 0.95, method = :woolf)
    @test collect(ci)  ≈ [0.2950420027, 0.9727082695] atol=1E-6

    # awoolf
    # Validate
    ci = MetidaFreq.orci(ct; level = 0.95, method = :awoolf)
    @test collect(ci)  ≈ [0.2982065561622649, 0.975836295197442] atol=1E-6

    # mover
    # Validate
    ci = MetidaFreq.orci(ct; level = 0.95, method = :mover)
    @test collect(ci)  ≈ [0.2963748435372293, 0.9689058534780502] atol=1E-6

end

@testset "  Risk ratio Confidence Intarvals                          " begin
    ct = MetidaFreq.contab([30 70; 40 50])

    # Miettinen & Nurminen
    # 0.4605492931 - 0.9820955908
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :mn)
    @test collect(ci)  ≈ [0.4605492931511954, 0.9820955908214944] atol=1E-6
    @test_nowarn MetidaFreq.rrci(0, 100, 1, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.rrci(1, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.rrci(0, 100, 0, 100; level = 0.95, method = :mn)
    @test_nowarn MetidaFreq.rrci(100, 100, 100, 100; level = 0.95, method = :mn)
    # FM
    # validation
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [0.46101548213819626, 0.981136152040161] atol=1E-6
    @test_nowarn MetidaFreq.rrci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.rrci(1, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.rrci(0, 100, 0, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.rrci(100, 100, 100, 100; level = 0.95, method = :fm)

    # cli # validation
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :cli)
    @test collect(ci)  ≈ [0.4663950370893218, 0.9860541079252757] atol=1E-6

    # li | Wald
    # 0.4624671992 - 0.9852050064
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :li)
    @test collect(ci)  ≈ [0.46246719923578444, 0.9852050064370166] atol=1E-6

    # mover
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :mover)
    @test collect(ci)  ≈ [0.46344425873893524, 0.9808806908109405] atol=1E-6

end

@testset "  Contab from tabular data                                 " begin
    freqdat  = CSV.File(path*"/csv/freqdat.csv") |> DataFrame

    ct = MetidaFreq.freq(freqdat, :row)
    @test ct.tab[1,1] == 26
    @test ct.tab[2,1] == 17
    @test size(ct) == (2,1)

    ct = MetidaFreq.freq(freqdat[!, :row])
    @test ct.tab[1,1] == 26
    @test ct.tab[2,1] == 17
    @test size(ct) == (2,1)

    ctd = MetidaFreq.freqdict(freqdat[!, :row])
    @test ctd["w"] == 26
    @test ctd["q"] == 17

    ct = MetidaFreq.contab(freqdat, :row, :col)
    @test ct.tab[1,1] == 21
    @test ct.tab[1,2] == 5
    @test ct.tab[2,1] == 8
    @test ct.tab[2,2] == 9
    @test size(ct) == (2,2)

    cts = MetidaFreq.sumrows(x-> x + 1, ct; coln = "Val")
    @test cts[1,1] == ct[1,1] + ct[1,2] + 2
    @test cts[2,1] == ct[2,1] + ct[2,2] + 2 

    ctac = MetidaFreq.addcol(ct, [34, 56]; coln = "NewCol")
    @test ctac[1,3] == 34
    @test ctac[2,3] == 56

    ctac = MetidaFreq.addcol(sum, ct; coln = "NewCol")
    @test ctac[1,3] == ct[1,1] + ct[1,2] 
    @test ctac[2,3] == ct[2,1] + ct[2,2]

    ctac = MetidaFreq.addcol((x,y) -> sum(x) + y, ct, [3, 6]; coln = "NewCol")
    @test ctac[1,3] == ct[1,1] + ct[1,2] + 3 
    @test ctac[2,3] == ct[2,1] + ct[2,2] + 6

    ctds = MetidaFreq.DataSet([deepcopy(ct), MetidaFreq.contab([1 2; 3 4], rownames = ["w", "q"], colnames = ["a", "b"])])
    ctac = MetidaFreq.colreduce(sum, ctds; coln = nothing)
    @test ctac[1,1] == 26
    @test ctac[1,2] == 3
    @test ctac[2,1] == 17
    @test ctac[2,2] == 7

    push!(ctds, MetidaFreq.contab([0 0; 0 0], rownames = ["w", "q"], colnames = ["a", "b"]))
    @test length(ctds) == 3
    MetidaFreq.dropzeros!(ctds)
    @test length(ctds) == 2

    @test_nowarn show(io, ctds)

    @test_nowarn permutedims(ct)
    @test_nowarn MetidaFreq.contab(ct, [1], 1:2)
    @test_nowarn MetidaFreq.contab(ct, 1, 1:2)

    Base.show(io, ct)

    HypothesisTests.ChisqTest(ct)

    HypothesisTests.MultinomialLRTest(ct)

    HypothesisTests.FisherExactTest(ct)

    transform!(freqdat, :row => categorical, renamecols=false)
    transform!(freqdat, :col => categorical, renamecols=false)

    ct = MetidaFreq.contab(freqdat, :row, :col)
    @test ct.tab[1,1] == 8
    @test ct.tab[1,2] == 9
    @test ct.tab[2,1] == 21
    @test ct.tab[2,2] == 5

    freqdat2  = CSV.File(path*"/csv/ft.csv") |> DataFrame

    ct = MetidaFreq.contab(freqdat2, :row, :col)
    @test ct.tab[1,1] == 20
    @test ct.tab[1,2] == 24
    @test ct.tab[2,1] == 83
    @test ct.tab[2,2] == 44

    ct = MetidaFreq.contab(freqdat2, :row, :col; sort = :s1)

    @test ct[1].tab[1,1] == 0
    @test ct[1].tab[1,2] == 2
    @test ct[1].tab[2,1] == 65
    @test ct[1].tab[2,2] == 34

    @test ct[2].tab[1,1] == 20
    @test ct[2].tab[1,2] == 22
    @test ct[2].tab[2,1] == 18
    @test ct[2].tab[2,2] == 10

    ct = MetidaFreq.contab(freqdat2, :row, :col; sort = [:s1, :s2])
    @test  length(ct) == 6

    @test ct[6].id[:s1] == "e"
    @test ct[6].id[:s2] == "i"

    @test ct[6].tab[1,1] == 6
    @test ct[6].tab[1,2] == 12
    @test ct[6].tab[2,1] == 18
    @test ct[6].tab[2,2] == 10

    transform!(freqdat2, :row => categorical, renamecols=false)
    transform!(freqdat2, :col => categorical, renamecols=false)
    transform!(freqdat2, :s1 => categorical, renamecols=false)
    transform!(freqdat2, :s2 => categorical, renamecols=false)

    ct = MetidaFreq.contab(freqdat2, :row, :col; sort = [:s1, :s2])

    @test ct[5].tab[1,1] == 12
    @test ct[5].tab[1,2] == 6
    @test ct[5].tab[2,1] == 10
    @test ct[5].tab[2,2] == 18

    ct = MetidaFreq.contab([1 2; 3 4])
    Base.show(io, ct)
    ct = MetidaFreq.contab([1, 2, 3, 4])
    Base.show(io, ct)
end

@testset "  Meta proportions                                         " begin
    # Make DataSet
    ############################################################################
    pf1 = MetidaFreq.contab([15 8; 5 14])
    pf2 = MetidaFreq.contab([45 72; 23 95])
    mds = MetidaFreq.DataSet([pf1, pf2])
    ############################################################################

    # Check all metrics
    mp = MetidaFreq.metaprop(mds, :rr)
    mp = MetidaFreq.metaprop(mds, :or)
    mp = MetidaFreq.metaprop(mds, :diff)
    show(io, mp)

    # Fixed effect MH (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    @test mpf.est ≈ 0.2196889005445779 atol=1E-6
    @test mpf.var ≈ 0.002912630389616186 atol=1E-6
    ci = confint(mpf; level = 0.95)
    @test collect(ci)  ≈ [0.11391201412030795, 0.32546578696884787] atol=1E-6
    show(io, mpf)

    # Fixed effect IV (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)

    # Check random effect tau (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metaproprandom(mp; tau = :dl)
    mpf = MetidaFreq.metaproprandom(mp; tau = :ho)
    mpf = MetidaFreq.metaproprandom(mp; tau = :hm)
    mpf = MetidaFreq.metaproprandom(mp; tau = :sj)
    show(io, mpf)


    pf1 = MetidaFreq.contab([30 70; 40 60])
    pf2 = MetidaFreq.contab([20 70; 30 60])
    pf3 = MetidaFreq.contab([30 80; 30 90])
    mds = MetidaFreq.DataSet([pf1, pf2, pf3])
    mp  = MetidaFreq.metaprop(mds, :diff)
    
    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)
    @test mpf.est ≈ -0.054583175 atol=1E-6
    @test sqrt(mpf.var) ≈ 0.036584324 atol=1E-6
    @test mpf.hetq ≈ 2.962180668 atol=1E-6
    ci = MetidaFreq.confint(mpf; level = 0.95)
    @test ci[1] ≈ -0.126287133 atol=1E-6
    @test ci[2] ≈  0.017120783 atol=1E-6

    #mpf.heti

    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)

    #
    mpf = MetidaFreq.metaproprandom(mp; tau = :dl)
    @test mpf.est ≈ -0.057406278 atol=1E-6
    @test sqrt(mpf.var) ≈ 0.044680842 atol=1E-6
    @test mpf.hetq ≈ 2.962180668 atol=1E-6
    @test mpf.heti ≈ 32.482173636 atol=1E-6
    @test mpf.hettau ≈ 0.001949933 atol=1E-6

    mpf = MetidaFreq.metaproprandom(mp; tau = :ho)
    @test mpf.est ≈ -0.056863933 atol=1E-6
    @test sqrt(mpf.var) ≈ 0.042684404 atol=1E-6
    @test mpf.hetq ≈ 2.962180668 atol=1E-6
    @test mpf.heti ≈ 26.097161948 atol=1E-6
    @test mpf.hettau ≈ 0.001431282 atol=1E-6

    mpf = MetidaFreq.metaproprandom(mp; tau = :hm)

    mpf = MetidaFreq.metaproprandom(mp; tau = :sj)
    #!!!!!
    #@test mpf.est ≈ -0.058065603 atol=1E-6
    #@test sqrt(mpf.var) ≈ 0.047558623 atol=1E-6
    #@test mpf.hetq ≈ 2.962180668 atol=1E-6
    #@test mpf.heti ≈ 40.340569062 atol=1E-6
    #@test mpf.hettau ≈ 0.002740665 atol=1E-6

    # 3 TRIAL CASE
# Validated against OpenMetaAnalyst
# http://www.cebm.brown.edu/openmeta/

    metadf  = CSV.File(path*"/csv/meta.csv") |> DataFrame
    ctds = MetidaFreq.contab(metadf, :group, :result; sort = :trial)

    # DIFF 
    mp = MetidaFreq.metaprop(ctds, :diff)
        # FIXED
        
#=
Binary Fixed-Effect Model - Inverse Variance
Metric: Risk Difference
 Model Results
 Estimate  Lower bound   Upper bound   Std. error   p-Value   
 0.141681   -0.011966      0.295327     0.078392    0.070712  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 4.423578    0.109505  
study names  weights   
trial 1: 23.428909%
trial 2: 34.782495%
trial 3: 41.788596%
=#
        mpf = MetidaFreq.metapropfixed(mp; weights = :iv)
        @test MetidaFreq.weights(mpf) ≈ [34.782495, 41.788596, 23.428909] atol=1E-5 
        ci = MetidaFreq.confint(mpf; level = 0.95)
        @test ci[1] ≈ -0.011966 atol=1E-5
        @test ci[2] ≈ 0.295327 atol=1E-5
        @test mpf.est ≈ 0.141681 atol=1E-5
        @test sqrt(mpf.var) ≈ 0.078392 atol=1E-5
        #@test mpf.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpf.hetq ≈ 4.423578 atol=1E-5

#=
Binary Fixed-Effect Model - Mantel Haenszel
Metric: Risk Difference
 Model Results
 Estimate  Lower bound   Upper bound   Std. error   p-Value   
 0.114863   -0.047325      0.277052     0.082751    0.165118  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 4.540603    0.103281 
 study names  weights   
trial 1: 27.010409%
trial 2: 26.572791%
trial 3: 46.416800%   
=#
#rma.mh(measure="RD", ai=c(21,14,11), bi=c(20,17,12), ci=c(2,12,9), di=c(10,22,7))

        mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
        @test MetidaFreq.weights(mpf) ≈ [26.572791, 46.416800, 27.010409] atol=1E-5 
        ci = MetidaFreq.confint(mpf; level = 0.95)
        @test ci[1] ≈ -0.047325 atol=1E-5
        @test ci[2] ≈ 0.277052 atol=1E-5
        @test mpf.est ≈ 0.114863 atol=1E-5
        @test sqrt(mpf.var) ≈ 0.082751 atol=1E-5
        #@test mpf.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpf.hetq ≈ 4.540603 atol=1E-5
        @test mpf.heti ≈ 55.95 atol=1E-2

        # RANDOM 
#=
Binary Random-Effects Model
Metric: Risk Difference
 Model Results
 Estimate  Lower bound   Upper bound   Std. error   p-Value   
 0.130758   -0.112013      0.373530     0.123865    0.291128  
 Heterogeneity
 tau^2     Q(df=2)    Het. p-Value      I^2     
 0.026981  4.423578     0.109505     58.776613 
 study names  weights   
trial 1: 28.833406%
trial 2: 34.362457%
trial 3: 36.804138%
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :ho)
        @test MetidaFreq.weights(mpr) ≈ [34.362457, 36.804138, 28.833406] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.112013 atol=1E-5
        @test ci[2] ≈ 0.373530 atol=1E-5
        @test mpr.est ≈ 0.130758 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.123865 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 4.423578 atol=1E-5
        @test mpr.heti ≈ 58.776613 atol=1E-2
        @test mpr.hettau ≈ 0.026981 atol=1E-2

#=
Binary Random-Effects Model
Metric: Risk Difference
 Model Results
 Estimate  Lower bound   Upper bound   Std. error   p-Value   
 0.131655   -0.100067      0.363378     0.118228    0.265464  
 Heterogeneity
 tau^2     Q(df=2)    Het. p-Value     I^2     
 0.022931  4.423578     0.109505     54.78773  
 study names  weights   
trial 1: 28.432718%
trial 2: 34.428819%
trial 3: 37.138462%
=#      
        mpr = MetidaFreq.metaproprandom(mp; tau = :dl)
        @test MetidaFreq.weights(mpr) ≈ [34.428819, 37.138462, 28.432718] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.100067 atol=1E-5
        @test ci[2] ≈ 0.363378 atol=1E-5
        @test mpr.est ≈ 0.131655 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.118228 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 4.423578 atol=1E-5
        @test mpr.heti ≈ 54.78773 atol=1E-2
        @test mpr.hettau ≈ 0.022931 atol=1E-2

#=
=#
        # NEED CHECK SOMEHOW
        mpr = MetidaFreq.metaproprandom(mp; tau = :hm)
        @test MetidaFreq.weights(mpr) ≈ [34.44537645338772, 37.22663989261683, 28.32798365399544] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.09715839439746077 atol=1E-5
        @test ci[2] ≈ 0.3609333913404578 atol=1E-5
        @test mpr.est ≈ 0.131887498471498535 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.11686229679506546 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 4.423577921113032 atol=1E-5
        @test mpr.heti ≈ 53.7359206749384 atol=1E-2
        @test mpr.hettau ≈ 0.021979683269097428 atol=1E-2

#=
Binary Random-Effects Model
Metric: Risk Difference
 Model Results
 Estimate  Lower bound   Upper bound   Std. error   p-Value   
 0.130641   -0.113672      0.374954     0.124652    0.294618  
 Heterogeneity
 tau^2     Q(df=2)    Het. p-Value      I^2     
 0.027562  4.423578     0.109505     59.291692  
 study names  weights   
trial 1: 28.885553%
trial 2: 34.353474%
trial 3: 36.760973%
=#

        mpr = MetidaFreq.metaproprandom(mp; tau = :sj)
#=
SPSS ML								
Effect size   SD		   Z		  Pval        95%CIL     95%CIU
0.136971	0.093606	1.463261	0.143396	-0.046495	0.320436
=#  
        # ML 
        mpr = MetidaFreq.metaproprandom(mp; tau = :ml)
        @test MetidaFreq.weights(mpr) ≈ [34.72637917278319, 39.345440640642835, 25.928180186573965] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.04649350144021322 atol=1E-5
        @test ci[2] ≈ 0.3204348488098385 atol=1E-5
        @test mpr.est ≈ 0.13697067368481264 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.09360589101236938 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 4.423577921113032 atol=1E-5
        @test mpr.heti ≈ 28.55612879933495 atol=1E-2
        @test mpr.hettau ≈ 0.007563712502569099 atol=1E-2

        # REML
        mpr = MetidaFreq.metaproprandom(mp; tau = :reml)
        @test MetidaFreq.weights(mpr) ≈ [34.424498831686016, 37.11579173325517, 28.459709435058826] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.10083053686105758 atol=1E-5
        @test ci[2] ≈ 0.364020576665633 atol=1E-5
        @test mpr.est ≈ 0.1315950199022877 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.11858664679386356 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 4.423577921113032 atol=1E-5
        @test mpr.heti ≈ 55.05817174207822 atol=1E-2
        @test mpr.hettau ≈ 0.023183110754144722 atol=1E-2

        
    # OR 

    mp = MetidaFreq.metaprop(ctds, :or)

        # FIXED 
#=
Binary Fixed-Effect Model - Inverse Variance
Metric: Odds Ratio
 Model Results
 Estimate  Lower bound   Upper bound   p-Value   
 1.516596    0.745928      3.083490    0.250011  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 3.540859    0.170260    
 Results (log scale)
 Estimate  Lower bound   Upper bound   Std. error  
 0.416468   -0.293126      1.126062     0.362044   
 study names  weights   
trial 1: 30.610197%
trial 2: 18.789070%
trial 3: 50.600734%
=#

    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)
    @test MetidaFreq.weights(mpf) ≈ [18.789070, 50.600734, 30.610197] atol=1E-5 
    ci = MetidaFreq.confint(mpf; level = 0.95)
    @test ci[1] ≈ -0.293126 atol=1E-5
    @test ci[2] ≈ 1.126062 atol=1E-5
    @test mpf.est ≈ 0.416468 atol=1E-5
    @test sqrt(mpf.var) ≈ 0.362044 atol=1E-5
    #@test mpf.chisq ≈ 4.864103956654892 atol=1E-5
    @test mpf.hetq ≈ 3.540859 atol=1E-5
    #@test mpf.heti ≈ 43.516535584304336 atol=1E-2 # CHECK
#=
Binary Fixed-Effect Model - Mantel Haenszel
Metric: Odds Ratio
 Model Results
 Estimate  Lower bound   Upper bound   p-Value   
 1.602286    0.812626      3.159290    0.173521  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 3.563907    0.168309    
 Results (log scale)
 Estimate  Lower bound   Upper bound   Std. error  
 0.471431   -0.207485      1.150347     0.346392
 study names  weights   
trial 1: 41.565005%
trial 2: 11.327989%
trial 3: 47.107006%
=#

    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    @test MetidaFreq.weights(mpf) ≈ [11.327989, 47.107006, 41.565005] atol=1E-5 
    ci = MetidaFreq.confint(mpf; level = 0.95)
    @test ci[1] ≈ -0.207485 atol=1E-5
    @test ci[2] ≈ 1.150347 atol=1E-5
    @test mpf.est ≈ 0.471431 atol=1E-5
    @test sqrt(mpf.var) ≈ 0.346392 atol=1E-5
    #@test mpf.chisq ≈  atol=1E-5
    @test mpf.hetq ≈ 3.563907 atol=1E-5
    #@test mpf.heti ≈ 43.88 atol=1E-2 # CHECK


        # RANDOM 

        mpr = MetidaFreq.metaproprandom(mp; tau = :ho)

#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :dl)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :hm)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :sj)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :ml)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :reml)

    # RR

    mp = MetidaFreq.metaprop(ctds, :rr)

        # FIXED
#=
Binary Fixed-Effect Model - Inverse Variance
Metric: Relative Risk
 Model Results
 Estimate  Lower bound   Upper bound   p-Value   
 1.161689    0.774870      1.741611    0.468192  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 3.266050    0.195338    
 Results (log scale)
 Estimate  Lower bound   Upper bound   Std. error  
 0.149875   -0.255060      0.554810     0.206603
 study names  weights   
trial 1: 44.444116%
trial 2:  9.703440%
trial 3: 45.852444%
=#
    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)
    @test MetidaFreq.weights(mpf) ≈ [9.703440, 45.852444, 44.444116] atol=1E-5 
    ci = MetidaFreq.confint(mpf; level = 0.95)
    @test ci[1] ≈ -0.255060 atol=1E-5
    @test ci[2] ≈ 0.554810 atol=1E-5
    @test mpf.est ≈ 0.149875 atol=1E-5
    @test sqrt(mpf.var) ≈ 0.206603 atol=1E-5
    #@test mpf.chisq ≈ 3.792288756192065 atol=1E-5
    @test mpf.hetq ≈ 3.266050 atol=1E-5
    #@test mpf.heti ≈ 38.76395140594738 atol=1E-2 # CHECK
#=
Binary Fixed-Effect Model - Mantel Haenszel
Metric: Relative Risk
 Model Results
 Estimate  Lower bound   Upper bound   p-Value   
 1.319025    0.869418      2.001139    0.192921  
 Heterogeneity
 Q(df=2)   Het. p-Value  
 3.644017    0.161701    
 Results (log scale)
 Estimate  Lower bound   Upper bound   Std. error  
 0.276893   -0.139931      0.693717     0.212669 
 study names  weights   
trial 1: 42.198426%
trial 2: 12.300662%
trial 3: 45.500912%
=#

    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    @test MetidaFreq.weights(mpf) ≈ [12.300662, 45.500912, 42.198426] atol=1E-5 
    ci = MetidaFreq.confint(mpf; level = 0.95)
    @test ci[1] ≈ -0.139931 atol=1E-5
    @test ci[2] ≈ 0.693717 atol=1E-5
    @test mpf.est ≈ 0.276893 atol=1E-5
    @test sqrt(mpf.var) ≈ 0.212669 atol=1E-5
    #@test mpf.chisq ≈ 3.792288756192065 atol=1E-5
    @test mpf.hetq ≈ 3.644017 atol=1E-5
    #@test mpf.heti ≈ 45.11551592477506 atol=1E-2 # CHECK

        # RANDOM 
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :ho)
        @test MetidaFreq.weights(mpr) ≈ [19.280363, 40.548276, 40.171359] atol=1E-5 
        ci = MetidaFreq.confint(mpr; level = 0.95)
        @test ci[1] ≈ -0.44854363729260716 atol=1E-5
        @test ci[2] ≈ 0.9510473341505784 atol=1E-5
        @test mpr.est ≈ 0.2512518484289856 atol=1E-5
        @test sqrt(mpr.var) ≈ 0.3570450739102805 atol=1E-5
        #@test mpr.chisq ≈ 7.690002992971811 atol=1E-5
        @test mpr.hetq ≈ 3.2660500569826847 atol=1E-5
        @test mpr.heti ≈ 60.172084286641514 atol=1E-2
        @test mpr.hettau ≈ 0.22130152193554511 atol=1E-2   
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :dl)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :hm)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :sj)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :ml)
#=
CHECK
=#
        mpr = MetidaFreq.metaproprandom(mp; tau = :reml)


    # DataFrame Test
    df = DataFrame(ctds)
    @test df[1, 1] == 2
    @test df[1, 2] == 21
    @test df[1, 3] == 20
    @test df[1, 4] == 2
    @test df[1, 5] == 10
end

@testset "  Goodman CI                                               " begin

    ci = MetidaFreq.ci_prop_goodman([44,55,43,32,67,78], 0.05)
    @test ci[1][1] ≈ 0.09468368035184335  atol=1E-6
    @test ci[1][2] ≈ 0.1966412812730604  atol=1E-6
    @test ci[6][1] ≈ 0.18692731624998904  atol=1E-6
    @test ci[6][2] ≈ 0.3130119423857655  atol=1E-6

    ci = MetidaFreq.ci_prop_goodman([91,49,37,43], 0.05)
    @test collect(ci[1]) ≈ [0.3342024783268433, 0.4978332073846343]  atol=1E-6
    @test collect(ci[3]) ≈ [0.11455134046157951, 0.24011208358778185]  atol=1E-6

    ct = MetidaFreq.contab(permutedims([91,49,37,43]))
    ci = MetidaFreq.mpropci(ct; level = 0.95, method = :default)
    @test collect(ci[1]) ≈ [0.3342024783268433, 0.4978332073846343]  atol=1E-6
    @test collect(ci[3]) ≈ [0.11455134046157951, 0.24011208358778185]  atol=1E-6

    ct = MetidaFreq.contab([91 49 37 43; 91 49 37 43])
    ci = MetidaFreq.mpropci(ct; level = 0.95)
    @test ci[1] == ci[2]
end


@testset "  Errors                                                   " begin
    ct = MetidaFreq.contab([38, 62])

    @test_throws ArgumentError MetidaFreq.propci(ct; level = 0.95, method = :err)
    @test_throws ArgumentError MetidaFreq.diffci(ct; level = 0.95, method = :err)
    @test_throws ArgumentError MetidaFreq.orci(ct; level = 0.95, method = :err)
    @test_throws ArgumentError MetidaFreq.rrci(ct; level = 0.95, method = :err)

    @test_throws ArgumentError MetidaFreq.diffci(ct; level = 0.95, method = :wald)
    @test_throws ArgumentError MetidaFreq.orci(ct; level = 0.95, method = :wald)
    @test_throws ArgumentError MetidaFreq.rrci(ct; level = 0.95, method = :wald)
end
