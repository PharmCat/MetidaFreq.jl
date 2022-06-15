using Documenter, MetidaFreq, Weave, PrettyTables, CSV, DataFrames
#using DocumenterLaTeX


makedocs(
        modules = [MetidaFreq],
        sitename = "MetidaFreq.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            "Contingency tables" => "contab.md",
            "Confidence intervals" => "ci.md",
            "HypothesisTests" => "ht.md",
            "Meta-analysis" => "meta.md",
            "Examples" => "examples.md",
            "Details" => "details.md",
            "API" => "api.md"
            ],
        )


deploydocs(repo = "github.com/PharmCat/MetidaFreq.jl.git", devbranch = "main", forcepush = true
)
