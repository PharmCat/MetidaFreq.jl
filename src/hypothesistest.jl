
function HypothesisTests.ChisqTest(contab::ConTab)
    HypothesisTests.ChisqTest(contab.tab)
end

function HypothesisTests.MultinomialLRTest(contab::ConTab)
    HypothesisTests.MultinomialLRTest(contab.tab)
end

function HypothesisTests.FisherExactTest(contab::ConTab)
    if !(size(contab, 1) == size(contab, 2) == 2)  error("FisherExactTest can be applied only for 2X2 tables") end
    HypothesisTests.FisherExactTest(contab.tab[1,1], contab.tab[1,2], contab.tab[2,1], contab.tab[2,2])
end
