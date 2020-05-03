using Documenter
using IntervalArithmetic, Turmeric, NumberIntervals

makedocs(
    modules = [Turmeric],
    doctest = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Turmeric.jl",
    authors = "David P. Sanders and Luis Benet, and Josua Grawitter",
    pages = Any[
        "Home" => "index.md",
        "`roots` interface" => "roots.md",
        "Internals" => "internals.md",
        "Bibliography" => "biblio.md",
        "API" => "api.md"
        ]
    )

deploydocs(
    repo = "github.com/gwater/Turmeric.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
