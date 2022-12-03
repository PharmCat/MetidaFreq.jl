# Examples

```@setup freqexample
using DataFrames
```

## Simple contingency table

```@example freqexample
using MetidaFreq, CSV, DataFrames;

df = CSV.File(joinpath(dirname(pathof(MetidaFreq)), "..", "test", "csv",  "ft.csv")) |> DataFrame

ct = MetidaFreq.contab(df, :row, :col)
```

## Confidence Intervals

```@example freqexample
MetidaFreq.propci(ct, method = :cp)
```

```@example freqexample
MetidaFreq.diffci(ct)
```

```@example freqexample
MetidaFreq.orci(ct)
```

```@example freqexample
MetidaFreq.rrci(ct)
```

## Sorting example

```@example freqexample
ct = MetidaFreq.contab(df, :row, :col; sort = :s1)
```

```@example freqexample
ct = MetidaFreq.contab(df, :row, :col; sort = [:s1, :s2])
```

```@example freqexample
MetidaFreq.dropzeros!(ct)
```

```@example freqexample
MetidaFreq.dropzeros!(ct)
```

## Meta-analysis

```@example freqexample
pf1 = MetidaFreq.contab([15 8; 5 14])
pf2 = MetidaFreq.contab([45 72; 23 95])
mds = MetidaFreq.DataSet([pf1, pf2])
```

```@example freqexample
mp = MetidaFreq.metaprop(mds, :rr)
```

```@example freqexample
mpf = MetidaFreq.metapropfixed(mp; weights = :mh)
```

```@example freqexample
mpf = MetidaFreq.metaproprandom(mp; tau = :dl)
```
