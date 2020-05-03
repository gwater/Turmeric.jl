
"""
    roots(f, X, method = Krawczyk(), tol = 1e-7, maxtasks = nothing)

Find all isolated roots of a function `f:R^n → R^n` in a region `X`, if the
number of roots is finite.

Inputs:
  - `f`: function whose roots will be found
  - `X`: `Interval` or `IntervalBox` in which roots are searched
  - `method`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol
    indicating the status. Current possible values are `Bisection()`, `Newton()`
    and `Krawczyk()`
  - `tol`: Absolute tolerance. If a region smaller than `tol` cannot be proven
    to contain or not contain a root, it is returned as indeterminate.
  - `maxtasks`: Limit the number of tasks spawned concurrently. By default set
    to `10^4 * nthreads()`.
"""
roots(f, X, method=Krawczyk(), tol=1e-7, maxtasks = nothing) =
    troots(f, X, method, tol, maxtasks)
