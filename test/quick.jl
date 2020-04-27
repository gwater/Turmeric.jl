using IntervalRootFinding, StaticArrays
include("../examples/smiley_examples.jl")
reg = SVector(SmileyExample52.region...)
IntervalRootFinding._roots(SmileyExample52.f, reg, Newton, DepthFirstSearch, 1e-5) |> display

using NumberIntervals
reg = NumberInterval.(SVector(SmileyExample52.region...))
IntervalRootFinding._roots(SmileyExample52.f, reg, Newton, DepthFirstSearch, 1e-5) |> display

