
export roots

function _roots(f, region, contractor, tol = default_tolerance, target_task_number = nothing)
    search = ThreadedRootSearch(f, region, contractor, tol, target_task_number)
    return last(search)
end

"""
    roots(f, region, method = Krawczyk(), tol = 1e-7, maxtasks = nothing)

Find all isolated roots of a function `f:R^n → R^n` in a region `region`, if the
number of roots is finite.

Inputs:
  - `f`: function whose roots will be found
  - `region`: region in which roots are searched given by an interval or a
    vector of intervals for multidimensional problems
  - `method`: function that, when applied to the function `f`, determines
    the status of a given box `region`. It returns the new box and a symbol
    indicating the status. Current possible values are `Bisection()`, `Newton()`
    and `Krawczyk()`
  - `tol`: Absolute tolerance. If a region smaller than `tol` cannot be proven
    to contain or not contain a root, it is returned as indeterminate.
  - `maxtasks`: Limit the number of tasks spawned concurrently. By default set
    to `10^4 * nthreads()`.
"""
function roots(
    f,
    region,
    method::Union{Krawczyk,Newton} = default_contractor,
    tol = default_tolerance,
    target_task_number = nothing
)
    contractor = GradientContractor(f, method, region)
    return _roots(f, region, contractor, tol, target_task_number)
end

function roots(
    f,
    region,
    method::Bisection,
    tol = default_tolerance,
    target_task_number = nothing
)
    return _roots(f, region, TrivialContractor(), tol, target_task_number)
end
