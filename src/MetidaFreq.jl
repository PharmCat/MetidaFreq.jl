# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

__precompile__()

module MetidaFreq

#using StaticArrays
#using AxisArrays
using Tables, CategoricalArrays, Distributions, Roots, StatsBase, LinearAlgebra, Optim

import HypothesisTests, MetidaBase
import StatsBase: confint
import MetidaBase: AbstractData, AbstractIdData, DataSet, Proportion, PrettyTables, metida_table_, getid, map, getdata
import HypothesisTests: ChisqTest, MultinomialLRTest, FisherExactTest
import Base: ht_keyindex, size, show, permutedims

export contab, propci, diffci, orci, rrci, confint

export metaprop, metapropfixed, metaproprandom

export ChisqTest, MultinomialLRTest, FisherExactTest

include("contab.jl")
include("confint.jl")
include("freq.jl")
include("metaprop.jl")
include("hypothesistest.jl")

end
