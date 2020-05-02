# IntervalRootFinding2.jl

[![Build Status](https://travis-ci.org/gwater/IntervalRootFinding2.jl.svg?branch=master)](https://travis-ci.org/gwater/IntervalRootFinding2.jl)

[![codecov.io](http://codecov.io/github/gwater/IntervalRootFinding2.jl/coverage.svg?branch=master)](http://codecov.io/github/gwater/IntervalRootFinding2.jl?branch=master)

This package provides guaranteed methods for finding **roots** of functions, i.e. solutions to the equation `f(x) == 0` for a function `f`.
To do so, it uses methods from interval analysis, using interval arithmetic from the [`IntervalArithmetic.jl`](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package by the same authors and ambiguity detection from [`NumberIntervals.jl`](https://github.com/gwater/NumberIntervals.jl) by Josua Grawitter.

## Basic usage examples

The basic function is `roots`. A standard Julia function and an interval is provided and the `roots` function return a list of intervals containing *all* roots of the function located in the starting interval.

```jl
julia> using NumberIntervals, IntervalRootFinding2

julia> f(x) = sin(x) - 0.1*x^2 + 1
f (generic function with 1 method)

julia> roots(f, NumberInterval(-10, 10))
4-element Array{Root{NumberInterval{Float64}},1}:
 Root([3.14959, 3.1496], :unique)
 Root([-4.42654, -4.42653], :unique)
 Root([-1.08205, -1.08204], :unique)
 Root([-3.10682, -3.10681], :unique)
```

The `:unique` status tell us, in addition, that each listed region contains exactly one root. The other possible status is `:unknown`, which corresponds to intervals that may contain zero, one, or more roots - no guarantee is provided for these intervals.

## Original Authors (IntervalRootFinding.jl)
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas,
Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders),
Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

### Acknowledgements ###

Financial support is acknowledged from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114, and DGAPA-UNAM PAPIIT grants IN-117214 and IN-117117. LB acknowledges support through a *Cátedra Moshinsky* (2013).
