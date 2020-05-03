# Internals

This section describes some of the internal mechanism of the package and several ways to use them to customize a search.

## Branch and bound

When `roots` is called, it performs a **branch-and-bound** search, that is, it iteratively looks at a region $X$ and for each region tries to determine if it contains a root. It then does the following:
  - If $X$ is proven to contain *no* root, it discards it.
  - If $X$ is proven to contain *exactly one root*, it tries to get the best possible bounds for the root and store the resulting region in the list of `Root`s to output with a `:unique` status.
  - If the test is inconclusive and the size of $X$ is smaller than the tolerance, it stores it in the list of `Root`s to output with `:unknown` status.
  - If the test is inconclusive and the size of $X$ is larger than the tolerance, it bisects $X$ and then processes each resulting half.

At some point all regions will either have a determined status or be smaller than the tolerance, and the algorithm will halt and return all stored roots.

## Contractors

To determine the status of a region, the algorithm uses so-called *contractors*. A `Contractor` is a callable object built from a function (in the case of `Bisection`) and possibly its derivative as well (for `Newton` and `Krawczyk`). When called with a region (wrapped in a `Root` object) and a tolerance, it returns the status of the root and the region (refined if the region contained a unique root).

```jl
julia> C = Newton(sin, cos)
Newton{typeof(sin),typeof(cos)}(sin, cos)

julia> C(Root(pi ± 0.001, :unknown), 1e-10)
Root([3.14159, 3.1416], :unique)

julia> C(Root(2 ± 0.001, :unkown), 1e-10)
Root([1.99899, 2.00101], :empty)
```

Contractors play a central role in the algorithm: they are the only part of it that varies for different methods.
