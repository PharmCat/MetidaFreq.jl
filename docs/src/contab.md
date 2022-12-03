# Contingency tables

Contingency tables can be constructed with [`contab`](@ref) method:

```
using MetidaFreq
contab(metadf, :group, :result; sort = :trial)
```

Where `metadf` - `DataFrame` or `Table`,  `:group` - row, `:result` - column. 

If you need separate tablel by value of some factor(s) - you can define `sort` value.