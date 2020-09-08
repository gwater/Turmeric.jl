
using NumberIntervals
using ThreadPools
using Turmeric
import Turmeric: GradientContractor, Newton, _isempty,
    strict_isinterior, bisect, contains_root, default_contractor, default_tolerance

import Base: push!, append!

struct StoreApplied{T, F}
    output::T
    f::F
end

push!(d::StoreApplied, x) = StoreApplied(push!(d.output, d.f(x)), d.f)
append!(d::StoreApplied, xs) = foreach(x -> push!(d, x), xs)

function sieve!(f, discarded)
    return function ff(data)
        append!(discarded, filter(f, data))
        return filter(!f, data)
    end
end

using Base.Threads

function mytmap(f, data)
    Ts = Base.return_types(f, (eltype(data),))
    out = Vector{first(Ts)}(undef, length(data))
    @threads for i in 1:length(data)
        out[i] = f(data[i])
    end
    return out
end

_tmap(f::Function) = data -> mytmap(f, data)
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

function setup_search(f, contractor, root_regions, indeterminate_regions, tolerance)
    with_contains_root = _appended(strictly_contains_root(f))
    with_contraction_and_isinterior =
        _appended(strictly_interior) ∘ _appended(contractor)
    with_region_too_small = _appended(too_small(tolerance) ∘ first)

    return function filter_contract_and_bisect(regions)
        remaining_region_tuples =
            regions |>
            _tmap(with_contains_root) |>
            _filter(last) |>
            _tmap(with_contraction_and_isinterior ∘ first) |>
            sieve!(last, StoreApplied(root_regions, first ∘ first)) |>
            _tmap(with_region_too_small ∘ first) |>
            sieve!(last, StoreApplied(indeterminate_regions, first ∘ first)) |>
            _tmap(_appended(!_isempty) ∘ _intersect ∘ first) |>
            _filter(last) |>
            _tmap(bisect ∘ first)
        return vcat(
            first.(remaining_region_tuples),
            last.(remaining_region_tuples)
        )
    end
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
    n_tasks = isnothing(num_tasks_hint) ? 4Threads.nthreads() :
        max(num_tasks_hint, Threads.nthreads())
    while length(regions) < n_tasks
        split_region_tuples = bisect.(regions)
        regions = vcat(
            first.(split_region_tuples),
            last.(split_region_tuples)
        )
    end
    # TODO use fold() method here
    while length(regions) > 0
        @show length(regions)
        regions = vcat(
            filter_contract_and_bisect(Iterators.take(regions, n_tasks) |> collect),
            @view regions[n_tasks + 1:end]
        )
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
include("../examples/smiley_examples.jl")
using .SmileyExample55

function main()
    region = NumberInterval(-50pi, 50pi)
    f(x) = sin(10sin(20sin(5x)))
    return roots(f, region)
end

function main2()
    region = NumberInterval.(SmileyExample55.region)
    return roots(SmileyExample55.f, region, Krawczyk(), 1e-6)
end

@time res = main2()[1]
@show length(res)
