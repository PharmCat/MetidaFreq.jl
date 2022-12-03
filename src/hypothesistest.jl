"""
    HypothesisTests.ChisqTest(contab::ConTab)
"""
function HypothesisTests.ChisqTest(contab::ConTab)
    HypothesisTests.ChisqTest(contab.tab)
end

"""
    HypothesisTests.MultinomialLRTest(contab::ConTab)
"""
function HypothesisTests.MultinomialLRTest(contab::ConTab)
    HypothesisTests.MultinomialLRTest(contab.tab)
end

"""
HypothesisTests.FisherExactTest(contab::ConTab)
"""
function HypothesisTests.FisherExactTest(contab::ConTab)
    if !(size(contab, 1) == size(contab, 2) == 2)  error("FisherExactTest can be applied only for 2X2 tables") end
    HypothesisTests.FisherExactTest(contab.tab[1,1], contab.tab[1,2], contab.tab[2,1], contab.tab[2,2])
end
