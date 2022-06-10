using Documenter, MetidaFreq, Weave, PrettyTables, CSV, DataFrames
#using DocumenterLaTeX


makedocs(
        modules = [MetidaFreq],
        sitename = "MetidaFreq.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            "Examples" => "examples.md",
            "Details" => "details.md",
            "Parameters" => "parameters.md",
            "API" => "api.md"
            ],
        )


deploydocs(repo = "github.com/PharmCat/MetidaFreq.jl.git", devbranch = "main", forcepush = true
)
