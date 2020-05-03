# `roots` interface

## Methods

Three root-finding methods currently available through the `roots` interface are the following:
  - `Newton`;
  - `Krawczyk` (default);
  - `Bisection`

Both the Newton and Krawczyk methods can determine if a root is unique in an interval, at the cost of requiring that the function is differentiable. The bisection method has no such requirement, but can never guarantee the existence or uniqueness of a root.

The method used is given as the third (optional) argument of the `roots` function:

```jl
julia> using NumberIntervals, Turmeric

julia> roots(log, NumberInterval(-2, 2), Newton())
(NumberInterval{Float64}[x ∈ [0.999592, 1.00026]], NumberInterval{Float64}[])

julia> roots(log, NumberInterval(-2, 2), Krawczyk())
(NumberInterval{Float64}[x ∈ [0.997849, 1.00145]], NumberInterval{Float64}[])

julia> roots(log, NumberInterval(-2, 2), Bisection())
(NumberInterval{Float64}[], NumberInterval{Float64}[x ∈ [0.999999, 1.00001]])
```

Note that as shown in the example, the `log` function does not complain about being given an interval going outside of its domain. While this may be surprising, this is the expected behavior and no root will ever be found outside the domain of a function.

In dimension greater than one, the function of interest must return a `SVector`, a type provided by the `StaticArrays` package, but otherwise works in the same way as in the 1D case.

```jl
julia> using StaticArrays

julia> function f( (x, y) )
           return SVector(sin(x), cos(y))
       end
f (generic function with 1 method)

julia> roots(f, SVector(NumberInterval(-3, 3), NumberInterval(-3, 3)))
(SArray{Tuple{2},NumberInterval{Float64},1,2}[[x ∈ [-6.12414e-06, 6.22146e-06], x ∈ [1.56843, 1.57243]], [x ∈ [-1.18531e-16, 1.20398e-16], x ∈ [-1.5708, -1.57079]]], SArray{Tuple{2},NumberInterval{Float64},1,2}[])
```

## Tolerance

An absolute tolerance for the search may be specified as the last argument of the `roots` function, the default being `1e-15`. Currently a method must first be provided in order to be able to choose the tolerance.

```jl
julia> g(x) = sin(exp(x))
g (generic function with 1 method)

julia> roots(g, NumberInterval(0, 2), Newton())
(NumberInterval{Float64}[x ∈ [1.14403, 1.1462], x ∈ [1.78904, 1.84259]], NumberInterval{Float64}[])

julia> roots(g, NumberInterval(0, 2), Newton(), 1e-2)
(NumberInterval{Float64}[x ∈ [1.14403, 1.1462], x ∈ [1.78904, 1.84259]], NumberInterval{Float64}[])
```

A lower tolerance may greatly reduce the computation time, at the cost of an increased number of regions being classified as indeterminate:

```jl
julia> using BenchmarkTools

julia> h(x) = cos(x) * sin(1 / x)
h (generic function with 1 method)

julia> @btime roots(h, NumberInterval(0.05, 1.0), Newton())
  221.925 μs (639 allocations: 42.47 KiB)
(NumberInterval{Float64}[x ∈ [0.318309, 0.318311], x ∈ [0.159147, 0.159161], x ∈ [0.0515382, 0.0531049], x ∈ [0.106089, 0.106113], x ∈ [0.0795764, 0.0795785], x ∈ [0.0636617, 0.0636622]], NumberInterval{Float64}[])

julia> @btime roots(h, NumberInterval(0.05, 1.0), Newton(), 1e-2)
  192.129 μs (549 allocations: 36.11 KiB)
(NumberInterval{Float64}[x ∈ [0.318309, 0.318311], x ∈ [0.159147, 0.159161], x ∈ [0.0515382, 0.0531049], x ∈ [0.106089, 0.106113]], NumberInterval{Float64}[x ∈ [0.0570253, 0.0641614], x ∈ [0.0785458, 0.0856819]])

julia> @btime roots(h, NumberInterval(0.05, 1.0), Newton(), 1e-1)
  106.830 μs (264 allocations: 16.83 KiB)
(NumberInterval{Float64}[], NumberInterval{Float64}[x ∈ [0.283803, 0.319372], x ∈ [0.05, 0.107542], x ∈ [0.107541, 0.165989]])
```

The last example shows a case where the tolerance was too large to be able to isolate the roots in distinct regions.

!!! warning

    For a root `x` of some function, if the absolute tolerance is smaller than `eps(x)` i.e. if `tol + x == x`, `roots` may never be able to converge to the required tolerance and the function may get stuck in an infinite loop.
