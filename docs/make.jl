using Documenter
using IntervalArithmetic, IntervalRootFinding2, NumberIntervals

makedocs(
    modules = [IntervalRootFinding2],
    doctest = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "IntervalRootFinding2.jl",
    authors = "David P. Sanders and Luis Benet",
    pages = Any[
        "Home" => "index.md",
        "`roots` interface" => "roots.md",
        "Internals" => "internals.md",
        "Bibliography" => "biblio.md",
        "API" => "api.md"
        ]
    )

deploydocs(
    repo = "github.com/gwater/IntervalRootFinding2.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
