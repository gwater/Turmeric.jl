using Base.Threads
using ThreadPools
using LazyArrays

import Base: last, iterate
export last, iterate

export ThreadedRootSearch

struct BisectionLimit end

struct ThreadBuffer{T, V <: AbstractVector{T}}
    root_regions::V
    indeterminate_regions::V
end
ThreadBuffer(region::T) where T = ThreadBuffer(T[], T[])

function refine!(buffer, region, contractor, f, tol)
    image = f(region)
    if !contains_root(image) || any(isempty.(image))
        return empty.(region)
    end

    contraction = contractor(region)
    if !any(isempty.(contraction)) # how would this happen?
        if isdisjoint(contraction, region)
            return empty.(region)
        end

        # UNIQUE?
        if strict_isinterior(contraction, region)
            push!(buffer.root_regions, contraction)
            return empty.(region)
        end

        region = region .∩ contraction # contract a bit
    else
        error()
    end

    if diam(region) < tol
        push!(buffer.indeterminate_regions, region)
        return empty.(region)
    end

    return region
end

function troots!(buffers, region, contractor, f, tol, generation, maxgenerations)

    buffer = buffers[threadid()]

    _region = refine!(buffer, region, contractor, f, tol)
    if any(isempty.(_region))
        return nothing
    end

    # RECURSION LIMIT?
    if generation > maxgenerations
        push!(buffer.indeterminate_regions, _region)
        return BisectionLimit()
    end

    # BISECT.
    a, b = bisect(_region)
    task1 = Threads.@spawn troots!(buffers, a, contractor, f, tol, generation + 1, maxgenerations)
    task2 = Threads.@spawn troots!(buffers, b, contractor, f, tol, generation + 1, maxgenerations)
    return task1, task2
end

isfinished(::Nothing) = true
isfinished(::BisectionLimit) = false
isfinished(t::Task) = isfinished(fetch(t))
isfinished(ts::Tuple{Task, Task}) = isfinished(ts[1]) & isfinished(ts[2])
isfinished(ts::AbstractVector) = mapreduce(isfinished, &, ts)

struct ThreadedRootSearch{B}
    buffers::B
    contractor
    f
    tol::Float64
    target_task_number::Int
    function ThreadedRootSearch(
        f,
        region::T,
        contractor,
        tol = default_tolerance,
        target_task_number = nothing
    ) where T
        n = nthreads()
        buffers = map(i -> ThreadBuffer(region), Tuple(1:n))
        push!(buffers[1].indeterminate_regions, region)
        return new{typeof(buffers)}(
            buffers,
            contractor,
            f,
            tol,
            ifelse(isnothing(target_task_number), 1_000n, target_task_number)
        )
    end
end

Base.IteratorSize(::Type{ThreadedRootSearch}) = Base.SizeUnknown()
Base.eltype(::Type{ThreadedRootSearch{B}}) where {N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}} =
    NTuple{2, ApplyVector{T, typeof(vcat), NTuple{N, V}}}

function chunk_regions!(buffers::NTuple{N}) where N
    vec_regions = map(buffers) do buffer
        unfinished_regions = copy(buffer.indeterminate_regions)
        # remove unfinished regions from buffers without freeing memory
        empty!(buffer.indeterminate_regions)
        return unfinished_regions
    end
    regions = reduce(vcat, vec_regions)
    n = length(regions)
    chunks = map(1:N:n) do i
        return @view regions[i:min(n, i + N - 1)]
    end
    return chunks
end

function iterate(
    search::ThreadedRootSearch{B},
    finished = false
) where {N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}}
    finished && return nothing
    chunks = chunk_regions!(search.buffers)
    finished = true
    for regions in chunks
        max_generations =
            max(1, floor(Int, log2(search.target_task_number / length(regions))))
        # trigger task cascade
        tasks = qmap(regions) do _region
            return troots!(
                search.buffers,
                _region,
                search.contractor,
                search.f,
                search.tol,
                1,
                max_generations
            )
        end
        # wait for cascade to finish
        finished &= isfinished(tasks)
    end
    root_regions = ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.root_regions, search.buffers)
    )
    indeterminate_regions =  ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.indeterminate_regions, search.buffers)
    )
    return (root_regions, indeterminate_regions), finished
end

function last(search::T) where T <: ThreadedRootSearch
    local res
    for s in search
        res = s
    end
    return res
end
