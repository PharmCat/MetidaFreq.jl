using MetidaFreq
using Test
using DataFrames, CSV, CategoricalArrays

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
    # SOC
    ci = MetidaFreq.propci(ct; level = 0.95, method = :soc)
    @test collect(ci)  ≈ [0.2891917018839231, 0.4775592393403469] atol=1E-6
    ci = MetidaFreq.propci(38, 100; level = 0.95, method = :soc)
    @test collect(ci)  ≈ [0.2891917018839231, 0.4775592393403469] atol=1E-6
    # Blaker - check!
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
    # Wald
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :wald)
    @test collect(ci)  ≈ [-0.28084842238, -0.00804046650] atol=1E-6
    # Wald CC
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :waldcc)
    @test collect(ci)  ≈ [-0.29140397794, 0.00251508905] atol=1E-6
    # Newcombe Hybrid Score
    ci = MetidaFreq.diffci(ct; level = 0.95, method = :nhs)
    @test collect(ci)  ≈ [-0.275381800, -0.007158419] atol=1E-6
    # # Newcombe Hybrid Score CC

end

@testset "  Odd ratio Confidence Intarvals                           " begin
end

@testset "  Risk ratio Confidence Intarvals                          " begin
end

@testset "  Contab from tabular data                                 " begin
    freqdat  = CSV.File(path*"/csv/freqdat.csv") |> DataFrame
    ct = MetidaFreq.contab(freqdat, :row, :col)
    @test ct.tab[1,1] == 21
    @test ct.tab[1,2] == 5
    @test ct.tab[2,1] == 8
    @test ct.tab[2,2] == 9

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
end
