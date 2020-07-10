
using NumberIntervals
using ThreadPools
using Turmeric
import Turmeric: GradientContractor, Newton, _isempty,
    strict_isinterior, bisect, contains_root, default_contractor, default_tolerance

import Base: push!, append!

struct StoreFirst{T}
    output::T
end

push!(d::StoreFirst, x) = StoreFirst(push!(d.output, first(x)))
append!(d::StoreFirst, xs) = foreach(x -> push!(d, x), xs)

function sieve!(f, discarded)
    return function ff(data)
        append!(discarded, filter(f, data))
        return filter(!f, data)
    end
end

_tmap(f::Function) = data -> tmap(f, data)
_filter(f::Function) = data -> filter(f, data)
_appended(f) = x -> (x, f(x))

strict_contains_zero(image) = contains_root(image) && !_isempty(image)
strictly_contains_root(f) = strict_contains_zero ∘ f
strictly_interior((region, contraction)) = strict_isinterior(contraction, region)
too_small(tol) = <(tol) ∘ diam
_intersect((region, contraction)::Tuple{T, T}) where T <: Number =
    region ∩ contraction
_intersect((region, contraction)::Tuple{T, T}) where T <: AbstractArray =
    region .∩ contraction

function _setup_search(with_contains_root, with_contraction_and_isinterior, with_too_small, root_regions, indeterminate_regions, tolerance)
    return function filter_contract_and_bisect(regions)
        remaining_region_tuples =
            regions |>
            _tmap(with_contains_root) |>
            _filter(last) |>
            _tmap(with_contraction_and_isinterior ∘ first) |>
            sieve!(last, StoreFirst(StoreFirst(root_regions))) |>
            _tmap(_appended(too_small(tolerance) ∘ first) ∘ first) |>
            sieve!(last, StoreFirst(StoreFirst(indeterminate_regions))) |>
            _tmap(bisect ∘ _intersect ∘ first)
        return vcat(
            first.(remaining_region_tuples),
            last.(remaining_region_tuples)
        )
    end
end

function setup_search(f, contractor, root_regions, indeterminate_regions, tolerance)
    with_contains_root = _appended(strictly_contains_root(f))
    with_contraction_and_isinterior =
        _appended(strictly_interior) ∘ _appended(contractor)
    with_too_small = _appended(too_small(tolerance) ∘ first)
    return _setup_search(with_contains_root, with_contraction_and_isinterior, with_too_small, root_regions, indeterminate_regions, tolerance)
end

function _roots(
    f,
    region,
    contractor,
    tolerance,
    num_tasks_hint
)
    indeterminate_regions = typeof(region)[]
    root_regions = typeof(region)[]
    regions = [region]
    filter_contract_and_bisect = setup_search(
        f,
        contractor,
        root_regions,
        indeterminate_regions,
        tolerance
    )
    n = isnothing(num_tasks_hint) ? 4Threads.nthreads() :
        max(num_tasks_hint, Threads.nthreads())
    while length(regions) > 0
        #@show length(regions)
        @show regions
        while length(regions) < n
            split_region_tuples = bisect.(regions)
            regions = vcat(
                first.(split_region_tuples),
                last.(split_region_tuples)
            )
        end
        regions = vcat(
            filter_contract_and_bisect(@view regions[1:n]),
            @view regions[n + 1:end]
        )
        sleep(2)
    end
    return root_regions, indeterminate_regions
end

function roots(
    f,
    region,
    method::Union{Krawczyk, Newton} = default_contractor,
    tolerance = default_tolerance,
    num_tasks_hint = nothing
)
    contractor = GradientContractor(f, method, region)
    return _roots(f, region, contractor, tolerance, num_tasks_hint)
end

# this method necessarily allocates memory in between steps
# we have multi threading in the map steps where there is immutability
# the sieve method mutates, so we run in single-threaded

# it might be more efficient to opportunistically calculate potentially
# unnecessary data to minimize sieve steps

# or we find a way to make sieve thread-safe

# not clear how we can avoid the allocations – we want to be non mutating

function main()
    f_lock = ReentrantLock()
    function f(x)
        lock(f_lock)
        res = sin(x)
        unlock(f_lock)
        return res
    end
    region = NumberInterval(-50pi, 50pi)
    return roots(f, region)
end

@time res = main()[1]
@show length(res)
