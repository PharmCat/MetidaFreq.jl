using MetidaFreq
using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();
freqdat  = CSV.File(path*"/csv/freqdat.csv") |> DataFrame



MetidaFreq.contab(freqdat, :row, :col)

freqtable(freqdat, :row, :col)

v = vcat(freqdat, freqdat)
for i = 1:10
    v = vcat(v, freqdat)
end


using Profile
using PProf
using BenchmarkTools


pprof()
PProf.kill()
Profile.clear()

@profile  for i = 1:10000  MetidaFreq.contab(v, :row, :col) end
@profile  for i = 1:10000  freqtable(v, :row, :col) end
@benchmark MetidaFreq.contab($v, :row, :col)
@benchmark freqtable($v, :row, :col)

freqdat = append!(freqdat, freqdat)
a = Vector{String}(undef, 86)
fill!(view(a, 1:40), "a")
fill!(view(a, 41:86), "b")
freqdat.s1 = a
freqdat = append!(freqdat, freqdat)
a = Vector{String}(undef, 172)
fill!(view(a, 1:100), "c")
fill!(view(a, 101:172), "d")
freqdat.s2 = a
MetidaFreq.contab(freqdat, :row, :col; sort = :s1)
MetidaFreq.contab(freqdat, :row, :col)

freqdat.s2 = Symbol.(freqdat.s2)
