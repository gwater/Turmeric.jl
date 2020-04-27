using IntervalRootFinding2, StaticArrays
include("../examples/smiley_examples.jl")
reg = SVector(SmileyExample52.region...)
IntervalRootFinding2._roots(SmileyExample52.f, reg, Newton, DepthFirstSearch, 1e-5) |> display

using NumberIntervals
reg = NumberInterval.(SVector(SmileyExample52.region...))
IntervalRootFinding2._roots(SmileyExample52.f, reg, Newton, DepthFirstSearch, 1e-5) |> display

