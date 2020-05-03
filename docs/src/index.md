# `Turmeric.jl`

This package provides guaranteed methods for finding **roots** of functions $f: \mathbb{R}^n \to \mathbb{R}^n$ with $n \ge 1$, i.e. vectors (or scalars, for $n=1$) $\mathbb{x}$ for which $f(\mathbb{x}) = \mathbb{0}$. In principle, it guarantees to find *all* roots inside a given box in $\mathbb{R}^n$, or report subboxes for which it is unable to provide guarantees.

To do so, it uses methods from interval analysis, using interval arithmetic from the [`IntervalArithmetic.jl`](https://github.com/JuliaIntervals/IntervalArithmetic.jl) and [`NumberIntervals.jl`](https://github.com/gwater/NumberIntervals.jl) packages.

!!! warning

    While this package aimed at providing *guaranteed results* and despite our best efforts and test suite, some bugs may remain and there are several known issues with corner cases. Please look at the [issue tracker](https://github.com/gwater/Turmeric.jl/issues) and report there any odd and/or incorrect behavior.

## Basic 1D example

To begin, we need a standard Julia function and an interval in which to search roots of that function. Intervals use the `NumberInterval` type provided by the `NumberIntervals.jl` package or the `Interval` type provided by `IntervalArithmetic.jl` and are generally constructed using the syntax, `NumberInterval(a, b)` (or `Interval(a, b)`) representing the closed interval $[a, b]$.

When provided with this information, the `roots` function will return a vector of all roots of the function in the given interval.

Example:

```jl
julia> using NumberIntervals, Turmeric

julia> rts = roots(x -> x^2 - 2x, NumberInterval(0, 10))
(NumberInterval{Float64}[x ∈ [1.97591, 2.0384]], NumberInterval{Float64}[x ∈ [0, 1.88237e-14]])
```

The roots are returned in two lists:
  - regions in the first list contain *exactly one* root of the function,
  - regions in the second list may or may not contain one or more roots; the algorithm used was unable to come to a conclusion.

The second status is still informative, since all regions of the original search interval *not* contained in *any* of the returned root intervals is guaranteed *not* to contain any root of the function. In the above example, we know that the function has no root in the interval $[2.1, 10]$, for example.

There are several known situations where the uniqueness (and existence) of a solution cannot be determined by the interval algorithms used in the package:
  - If the solution is on the boundary of the interval (as in the previous example);
  - If the derivative of the solution is zero at the solution.

In particular, the second condition means that multiple roots cannot be proven to be unique. For example:

```jl
julia> g(x) = (x^2 - 2)^2 * (x^2 - 3)
g (generic function with 1 method)

julia> roots(g, NumberInterval(-10, 10))
(NumberInterval{Float64}[x ∈ [1.70588, 1.75782], x ∈ [-1.73348, -1.73074]], NumberInterval{Float64}[x ∈ [-1.41422, -1.41421], x ∈ [1.41421, 1.41422]])
```

Here we see that the two double roots are reported as being possible roots without guarantee and the simple roots have been proved to be unique.


## Basic multi-dimensional example

For dimensions $n > 1$, the function passed to `roots` must currently return an `SVector` from the `StaticArrays.jl` package.

Here we give a 3D example:

```jl
julia> using StaticArrays

julia> function g( (x1, x2, x3) )
           return SVector(x1^2 + x2^2 + x3^2 - 1,
                          x1^2 + x3^2 - 0.25,
                          x1^2 + x2^2 - 4x3
                         )
       end
g (generic function with 1 method)

julia> X = NumberInterval(-5, 5)
x ∈ [-5, 5]

julia> rts = roots(g, SVector(X, X, X))[1]
4-element LazyArrays.ApplyArray{SArray{Tuple{3},NumberInterval{Float64},1,3},1,typeof(vcat),Tuple{Array{SArray{Tuple{3},NumberInterval{Float64},1,3},1}}}:
 [x ∈ [-0.441208, -0.440316], x ∈ [-0.866481, -0.865619], x ∈ [0.2359, 0.236236]]
 [x ∈ [-0.440977, -0.440598], x ∈ [0.865824, 0.866259], x ∈ [0.235966, 0.23617]]
 [x ∈ [0.44076, 0.440766], x ∈ [-0.866029, -0.866022], x ∈ [0.236066, 0.23607]]
 [x ∈ [0.440201, 0.441551], x ∈ [0.865395, 0.866755], x ∈ [0.235797, 0.236338]]
```

Thus the system admits four unique roots in the box $[-5, 5]^3$.
