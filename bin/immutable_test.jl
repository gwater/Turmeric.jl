
using NumberIntervals
using Turmeric
import Turmeric: GradientContractor, Newton, _isempty,
    strict_isinterior, bisect, contains_root

contains_roots((region, image)) = !_isempty(image) && contains_root(image)
strictly_interior((region, contraction)) = strict_isinterior(contraction, region)
append_contraction(contractor) = region -> (region, contractor(region))
region_too_small(tol) = region_contraction -> diam(first(region_contraction)) < tol
_intersect((region, contraction)) = region .∩ contraction
_appended(f) = x -> (x, f(x))

import ThreadPools: tmap
import Base: map

function sieve!(f, discarded)
    return function ff(data)
        append!(discarded, filter(f, data))
        return filter(!f, data)
    end
end
tmap(f::Function) = data -> tmap(f, data)
map(f::Function) = data -> map(f, data)

const _tmap = map

# this method necessarily allocates memory in between steps
# we have multi threading in the map steps where there is immutability
# the sieve method mutates, so we run in single-threaded

# it might be more efficient to opportunistically calculate potentially
# unnecessary data to minimize sieve steps

# or we find a way to make sieve thread-safe

# not clear how we can avoid the allocations – we want to be non mutating

import Base: push!, append!

struct StoreFirst{T}
    output::T
end

push!(d::StoreFirst, x) = StoreFirst(push!(d.output, first(x)))
append!(d::StoreFirst, xs) = foreach(x -> push!(d, x), xs)

function search!(f, contractor, root_regions, indeterminate_regions, regions, tol, n)
    remaining_region_tuples =
        Iterators.take(regions, n) |>
        _tmap(_appended(contains_roots) ∘ _appended(f)) |>
        data -> filter(last, data) |>
        _tmap(_appended(strictly_interior) ∘ append_contraction(contractor) ∘ first ∘ first) |>
        sieve!(last, StoreFirst(StoreFirst(root_regions))) |>
        _tmap(_appended(region_too_small(tol)) ∘ first) |>
        sieve!(last, StoreFirst(StoreFirst(indeterminate_regions))) |>
        _tmap(bisect ∘ _intersect ∘ first)
    remaining_regions = vcat(
        first.(remaining_region_tuples),
        last.(remaining_region_tuples),
        (@view regions[n+1:end])
    )
    return remaining_regions
end

function main()
    f = sin ∘ inv
    region = NumberInterval(-0.5pi, 0.5pi)
    indeterminate_regions = typeof(region)[]
    root_regions = typeof(region)[]
    regions = [region]
    while length(regions) > 0
        @show length(regions)
        regions = search!(
            f,
            GradientContractor(sin ∘ inv, Newton(), region),
            root_regions,
            indeterminate_regions,
            regions,
            1e-2,
            4Threads.nthreads()
        )
    end
    return root_regions, indeterminate_regions
end

@show main()[1]
