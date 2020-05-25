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

function refine_region!(buffer, region, contractor, f, tol)
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

refine_region!(buffer, region, search) =
    refine_region!(buffer, region, search.contractor, search.f, search.tol)

function complete!(search, region, generation, maxgeneration)
    buffer = search.buffers[threadid()]

    _region = refine_region!(buffer, region, search)
    if any(isempty.(_region))
        return true
    end

    if generation == maxgeneration
        append!(buffer.indeterminate_regions, bisect(_region))
        return false
    end

    tasks = map(bisect(_region)) do i
        Threads.@spawn complete!(search, i, generation + 1, maxgeneration)
    end
    # we cannot wait for results here because we need to free up the thread
    return tasks
end

complete!(b::Bool) = b
complete!(t::Task) = complete!(fetch(t))
complete!(ts::NTuple{N, Task}) where N = mapreduce(complete!, &, ts)
complete!(ts::AbstractVector) = mapreduce(complete!, &, ts)

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
    completed = false
) where {N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}}
    completed && return nothing
    chunks = chunk_regions!(search.buffers)
    completed = true
    for regions in chunks
        max_generations =
            max(1, floor(Int, log2(search.target_task_number / length(regions))))
        # trigger task cascade
        tasks = qmap(regions) do _region
            return complete!(search, _region, 1, max_generations)
        end
        # wait for cascade to finish
        completed &= complete!(tasks)
    end
    root_regions = ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.root_regions, search.buffers)
    )
    indeterminate_regions =  ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.indeterminate_regions, search.buffers)
    )
    return (root_regions, indeterminate_regions), completed
end

function last(search::T) where T <: ThreadedRootSearch
    local res
    for s in search
        res = s
    end
    return res
end
