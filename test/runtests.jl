using MetidaFreq
using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();
freqdat  = CSV.File(path*"/csv/freqdat.csv") |> DataFrame

#dat.li2007
#dat.hine1989
#dat.graves2010
#dat.bourassa1996
categorical!(freqdat, :row);
categorical!(freqdat, :col);

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
