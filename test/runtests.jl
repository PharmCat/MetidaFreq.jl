using MetidaFreq
using Test
using DataFrames, CSV, CategoricalArrays

path     = dirname(@__FILE__)
io       = IOBuffer();
freqdat  = CSV.File(path*"/csv/freqdat.csv") |> DataFrame

#dat.li2007
#dat.hine1989
#dat.graves2010
#dat.bourassa1996

transform!(freqdat, :row => categorical, renamecols=false)
transform!(freqdat, :col => categorical, renamecols=false)

@testset "MetidaFreq.jl" begin

    MetidaFreq.contab(freqdat, :row, :col)

    pf1 = MetidaFreq.contab([15 8; 5 14])
    pf2 = MetidaFreq.contab([45 72; 23 95])
    mds = MetidaFreq.DataSet([pf1, pf2])
    mp = MetidaFreq.metaprop(mds, :rr)
    mp = MetidaFreq.metaprop(mds, :or)
    mp = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    mpf = MetidaFreq.metapropfixed(mp; weights = :iv)
    mpf = MetidaFreq.metaproprandom(mp; tau = :dl)
    mpf = MetidaFreq.metaproprandom(mp; tau = :ho)
    mpf = MetidaFreq.metaproprandom(mp; tau = :hm)
    mpf = MetidaFreq.metaproprandom(mp; tau = :sj)

    pf1 = MetidaFreq.contab([15 8; 5 14])
    pf2 = MetidaFreq.contab([45 72; 23 95])
    mds = MetidaFreq.DataSet([pf1, pf2])
    mp = MetidaFreq.metaprop(mds, :diff)
    mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
    mpf = MetidaFreq.metaproprandom(mp; tau = :ho)

end

@testset "  Goodman CI" begin

    ci = MetidaFreq.ci_prop_goodman([44,55,43,32,67,78], 0.05)
    @test ci[1][1] ≈ 0.09468368035184335  atol=1E-6
    @test ci[1][2] ≈ 0.1966412812730604  atol=1E-6
    @test ci[6][1] ≈ 0.18692731624998904  atol=1E-6
    @test ci[6][2] ≈ 0.3130119423857655  atol=1E-6

    ci = MetidaFreq.ci_prop_goodman([91,49,37,43], 0.05)
    @test collect(ci[1]) ≈ [0.3342024783268433, 0.4978332073846343]  atol=1E-6
    @test collect(ci[3]) ≈ [0.11455134046157951, 0.24011208358778185]  atol=1E-6
end
