# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

__precompile__()

module MetidaFreq

#using StaticArrays
#using AxisArrays
using Tables, CategoricalArrays

import MetidaBase: AbstractData, DataSet

import Base.ht_keyindex

export contab

include("contab.jl")

end
