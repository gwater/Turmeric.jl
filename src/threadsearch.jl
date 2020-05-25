using Base.Threads
using ThreadPools
using LazyArrays

import Base: last, iterate
export last, iterate

export ThreadedRootSearch

struct ThreadBuffer{T, V <: AbstractVector{T}}
    root_regions::V
    indeterminate_regions::V
end
ThreadBuffer(region::T) where T = ThreadBuffer(T[], T[])

"""
    explore_region!(buffer, region, contractor, f, tol)

Explore `region` for roots of `f` using `contractor` and store findings in
`buffer`, up to search tolerance `tol`.
Returns parts of `region` which are left unexplored.
"""
function explore_region!(buffer, region, contractor, f, tol)
    image = f(region)
    if !contains_root(image) || any(isempty.(image))
        return empty.(region)
    end

    contraction = contractor(region)
    if strict_isinterior(contraction, region)
        push!(buffer.root_regions, contraction)
        return empty.(region)
    end

    if diam(region) < tol
        push!(buffer.indeterminate_regions, region)
        return empty.(region)
    end

    return region .∩ contraction
end

struct ThreadedRootSearch{B}
    buffers::B
    contractor
    f
    tol::Float64
    target_task_number::Int
    function ThreadedRootSearch(
        f,
        region,
        contractor,
        tol = default_tolerance,
        target_task_number = nothing,
        n = nothing
    )
        n = isnothing(n) ? nthreads() : n
        buffers = map(i -> ThreadBuffer(region), Tuple(1:n))
        push!(buffers[1].indeterminate_regions, region)
        return new{typeof(buffers)}(
            buffers,
            contractor,
            f,
            tol,
            isnothing(target_task_number) ? 1000n : target_task_number
        )
    end
end

Base.IteratorSize(::Type{ThreadedRootSearch}) = Base.SizeUnknown()
Base.eltype(::Type{ThreadedRootSearch{B}}) where {
    N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}
} = NTuple{2, ApplyVector{T, typeof(vcat), NTuple{N, V}}}

"""
    (search::ThreadedRootSearch)(region, generation, maxgeneration)

Advances the `search` by exploring and, if necessary, bisecting `region`.
Returns `true` when a region is explored up to search tolerance.
Returns `false` when the next step (`generation`) would exceed `maxgeneration`.
If exploration requires further bisection, returns an iterator of `Task`s.
Use `guarded_reduce(&, tasks)` to recursively fetch the boolean result.
"""
function (search::ThreadedRootSearch)(region, generation, maxgeneration)
    buffer = search.buffers[threadid()]

    _region = explore_region!(buffer, region, search.contractor, search.f, search.tol)
    if any(isempty.(_region))
        return true
    end

    if generation == maxgeneration
        append!(buffer.indeterminate_regions, bisect(_region))
        return false
    end

    tasks = map(bisect(_region)) do i
        Threads.@spawn search(i, generation + 1, maxgeneration)
    end
    # we cannot wait for results here because we need to free up the thread
    return tasks
end

guarded_reduce(op, b) = b
guarded_reduce(op, t::Task) = guarded_reduce(op, fetch(t))
guarded_reduce(op, ts::NTuple{N}) where N =
    mapreduce(t -> guarded_reduce(op, t), op, ts)
guarded_reduce(op, ts::AbstractVector) =
    mapreduce(t -> guarded_reduce(op, t), op, ts)

function pop_unfinished!(buffers)
    return mapreduce(vcat, buffers) do buffer
        unfinished_regions = copy(buffer.indeterminate_regions)
        # remove unfinished regions from buffers without freeing memory
        empty!(buffer.indeterminate_regions)
        return unfinished_regions
    end
end

function chunks(regions, chunksize)
    n = length(regions)
    return map(1:chunksize:n) do i
        return @view regions[i:min(n, i + chunksize - 1)]
    end
end

function collect_regions(buffers::B)  where {
    N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}
}
    root_regions = ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.root_regions, buffers)
    )
    indeterminate_regions =  ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.indeterminate_regions, buffers)
    )
    return root_regions, indeterminate_regions
end

estimate_generations(target_task_number, nregions) =
    max(1, floor(Int, log2(target_task_number / nregions)))

function iterate(search::ThreadedRootSearch{B}, completed = false) where {N, B <: NTuple{N}}
    completed && return nothing

    unfinished_regions = pop_unfinished!(search.buffers)
    completed = mapreduce(&, chunks(unfinished_regions, N)) do regions
        max_generations = estimate_generations(search.target_task_number, length(regions))
        tasks = qmap(r -> search(r, 1, max_generations), regions)
        return guarded_reduce(&, tasks)
    end

    return collect_regions(search.buffers), completed
end

function last(search::T) where T <: ThreadedRootSearch
    local res
    for s in search
        res = s
    end
    return res
end
