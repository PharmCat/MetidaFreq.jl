
using Test
using DataFrames, CSV, CategoricalArrays, HypothesisTests

path     = dirname(@__FILE__)
io       = IOBuffer();


#dat.li2007
#dat.hine1989
#dat.graves2010
#dat.bourassa1996


@testset "  Proportion Confidence Intarvals                          " begin
    ct = MetidaFreq.contab([38, 62])

    # Wald
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [0.28486600512143206, 0.47513399487856794] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [0.28486600512143206, 0.47513399487856794] atol=1E-6

    # Wald cc
    ci = MetidaFreq.propci(ct; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [0.27986600512143206, 0.48013399487856795] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [0.27986600512143206, 0.48013399487856795] atol=1E-6

    # Wilson
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wilson)
    @test collect(ci)  ≈ [0.2909759925247873, 0.47790244704488943] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wilson)
    @test collect(ci)  ≈ [0.2909759925247873, 0.47790244704488943] atol=1E-6

    # Wilson cc - check
    ci = MetidaFreq.propci(ct; level = 0.95, method = :wilsoncc)
    @test collect(ci)  ≈ [0.28639471634406627, 0.48294114679105093] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :wilsoncc)
    @test collect(ci)  ≈ [0.28639471634406627, 0.48294114679105093] atol=1E-6

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

    # FM | Mee
    # Validate
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [-0.27778778468650384, -0.007071207814437397] atol=1E-6
    @test_nowarn MetidaFreq.diffci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.diffci(1, 100, 0, 100; level = 0.95, method = :fm)

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
    # FM
    # 0.2958336891 - 0.9701570484
    ci = MetidaFreq.orci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [0.2958336891, 0.9701570484] atol=1E-6
    @test_nowarn MetidaFreq.orci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.orci(1, 100, 0, 100; level = 0.95, method = :fm)
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
    # FM
    # validation
    ci = MetidaFreq.rrci(ct; level = 0.95, method = :fm)
    @test collect(ci)  ≈ [0.46101548213819626, 0.981136152040161] atol=1E-6
    @test_nowarn MetidaFreq.rrci(0, 100, 1, 100; level = 0.95, method = :fm)
    @test_nowarn MetidaFreq.rrci(1, 100, 0, 100; level = 0.95, method = :fm)

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
    ct = MetidaFreq.contab(freqdat, :row, :col)
    @test ct.tab[1,1] == 21
    @test ct.tab[1,2] == 5
    @test ct.tab[2,1] == 8
    @test ct.tab[2,2] == 9

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

    # Fixed effect MH (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    @test mpf.est ≈ 0.2196889005445779 atol=1E-6
    @test mpf.var ≈ 0.0028728509817330943 atol=1E-6
    ci = confint(mpf; level = 0.95)
    @test collect(ci)  ≈ [0.1146368241999458, 0.32474097688921] atol=1E-6

    # Fixed effect IV (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)

    # Check random effect tau (diff)
    mp  = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metaproprandom(mp; tau = :dl)
    mpf = MetidaFreq.metaproprandom(mp; tau = :ho)
    mpf = MetidaFreq.metaproprandom(mp; tau = :hm)
    mpf = MetidaFreq.metaproprandom(mp; tau = :sj)

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
