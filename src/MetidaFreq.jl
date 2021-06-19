# Metida
# Copyright © 2019-2020 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>

__precompile__()

module MetidaFreq

#using StaticArrays
#using AxisArrays
using Tables, CategoricalArrays, Distributions, Roots, StatsBase

import HypothesisTests

import MetidaBase: AbstractData, AbstractIdData, DataSet, Proportion
import HypothesisTests: ChisqTest, MultinomialLRTest, FisherExactTest
import Base: ht_keyindex, size

export contab, diffci, confint

export metaprop, metapropfixed, metaproprandom

export ChisqTest, MultinomialLRTest, FisherExactTest

include("contab.jl")
include("confint.jl")
include("metaprop.jl")
include("hypothesistest.jl")

end
