# MetidaFreq.jl

This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

| Status | Cover | Build | Docs |
|--------|-------|-------|------|
|[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)|[![codecov](https://codecov.io/gh/PharmCat/MetidaFreq.jl/branch/main/graph/badge.svg?token=5Y2OJ9SJIE)](https://codecov.io/gh/PharmCat/MetidaFreq.jl)|[![Tier 1](https://github.com/PharmCat/MetidaFreq.jl/actions/workflows/Tier1.yml/badge.svg?branch=main)](https://github.com/PharmCat/MetidaFreq.jl/actions/workflows/Tier1.yml)|[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/MetidaFreq.jl/dev/) [![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://pharmcat.github.io/MetidaFreq.jl/stable/)|

Metida frequency tables.


## Installation

```
import Pkg; Pkg.add(url = "https://github.com/PharmCat/MetidaFreq.jl.git")
```

## Using

Contingency table from tabular data

```
MetidaFreq.contab(freqdat, :row, :col)
```

## See also

* [https://github.com/nalimilan/FreqTables.jl](https://github.com/nalimilan/FreqTables.jl)
