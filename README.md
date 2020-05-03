# IntervalRootFinding2.jl

[![Build Status](https://travis-ci.org/gwater/IntervalRootFinding2.jl.svg?branch=master)](https://travis-ci.org/gwater/IntervalRootFinding2.jl)

[![codecov.io](http://codecov.io/github/gwater/IntervalRootFinding2.jl/coverage.svg?branch=master)](http://codecov.io/github/gwater/IntervalRootFinding2.jl?branch=master)

This package provides guaranteed methods for finding **roots** of functions, i.e. solutions to the equation `f(x) == 0` for a function `f` using the multi-threading features first introduced in Julia 1.3.
To do so, it uses methods from interval analysis, using interval arithmetic from the [`IntervalArithmetic.jl`](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package and ambiguity detection from [`NumberIntervals.jl`](https://github.com/gwater/NumberIntervals.jl).

## Basic usage examples

The basic function is `roots`. A standard Julia function and an interval is provided and the `roots` function return a list of intervals containing *all* roots of the function located in the starting interval.

```jl
julia> using NumberIntervals, IntervalRootFinding2

julia> f(x) = sin(x) - 0.1*x^2 + 1
f (generic function with 1 method)

julia> roots(f, NumberInterval(-10, 10))[1]
4-element LazyArrays.ApplyArray{NumberInterval{Float64},1,typeof(vcat),Tuple{Array{NumberInterval{Float64},1}}}:
  x ∈ [3.07363, 3.25133]
 x ∈ [-1.14528, -1.01158]
 x ∈ [-4.43535, -4.41877]
 x ∈ [-3.10817, -3.10529]
```

The `roots()` function returns a tuple of two lists: The first lists contains all intervals which are proven to contain exactly one root of `f` and the second list contains other intervals which could not be excluded from the solution set.

The multi-threaded implementation was written by Josua Grawitter based on the
serial implementation in [IntervalRootFinding.jl](https://github.com/JuliaIntervals/IntervalRootFinding.jl).
