using Base.Threads
using ThreadPools
using LazyArrays

import Base: last
export last

export troots, tfullcontract, fullcontract, ThreadedRootSearch

struct BisectionLimit end

struct ThreadBuffer{T, V <: AbstractVector{T}}
    root_regions::V
    indeterminate_regions::V
end
ThreadBuffer(region::T) where T = ThreadBuffer(T[], T[])

function troots!(buffers, region, contractor, f, tol, generation, maxgenerations)

    buffer = buffers[threadid()]

    # DISCARD?
    image = f(region)
    if !contains_root(image) || isempty(image)
        return nothing
    end
    contraction = contractor(region)
    contraction_empty = any(isempty.(contraction))
    if isdisjoint(contraction, region) && !contraction_empty
        return nothing
    end

    # UNIQUE?
    if isinterior(contraction, region) && !contraction_empty
        push!(buffer.root_regions, contraction)
        return nothing
    end

    # TOLERANCE LIMIT?
    if diam(region) < tol
        push!(buffer.indeterminate_regions, region)
        return nothing
    end

    # CONTRACTABLE?
    if !contraction_empty
        region = region .∩ contraction
    end

    # RECURSION LIMIT?
    if generation > maxgenerations || Sys.free_memory() / 2^20 < 100
        push!(buffer.indeterminate_regions, region)
        return BisectionLimit()
    end

    # BISECT.
    a, b = bisect(region)
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
    maxgenerations::Int
    function ThreadedRootSearch(region::T, contractor, f, tol = 1e-7, maxgenerations = 20) where T
        n = nthreads()
        buffers = map(i -> ThreadBuffer(region), Tuple(1:n))
        push!(buffers[1].indeterminate_regions, region)
        return new{typeof(buffers)}(buffers, contractor, f, tol, maxgenerations)
    end
end

Base.IteratorSize(::Type{ThreadedRootSearch}) = Base.SizeUnknown()
Base.eltype(::Type{ThreadedRootSearch{B}}) where {N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}} =
    NTuple{2, ApplyVector{T, typeof(vcat), NTuple{N, V}}}

function iterate(
    search::ThreadedRootSearch{B},
    finished = false
) where {N, T, V, B <: NTuple{N, ThreadBuffer{T, V}}}
    finished && return nothing
    vec_regions = map(search.buffers) do buffer
        #crit = region -> diam(_region) > tol
        crit = region -> true
        unfinished_regions =
            filter(crit, buffer.indeterminate_regions)
        # remove unfinished regions from buffers
        filter!(!crit, buffer.indeterminate_regions)
        return unfinished_regions
    end
    regions = reduce(vcat, vec_regions)
    # trigger task cascade
    tasks = map(regions) do _region
        return troots!(
            search.buffers,
            _region,
            search.contractor,
            search.f,
            search.tol,
            1,
            search.maxgenerations
        )
    end
    root_regions = ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.root_regions, search.buffers)
    )
    indeterminate_regions =  ApplyVector{T, typeof(vcat), NTuple{N, V}}(
        vcat,
        map(b -> b.indeterminate_regions, search.buffers)
    )
    return (root_regions, indeterminate_regions), isfinished(tasks)::Bool
end

function last(search::T) where T <: ThreadedRootSearch
    local res
    for s in search
        res = s
    end
    return res
end

function troots(region, contractor, f, tol = 1e-7, maxgenerations = 20)
    search = ThreadedRootSearch(region, contractor, f, tol, maxgenerations)
    return last(search)
end

function fullcontract(region::T, contractor, tol = 1e-7, maxiters = 100) where T
    for i in 1:maxiters
        region2 = region .∩ contractor(region)
        if maximum(diam.(region2)) < maximum(diam.(region))
            region = region2
        else
            break
        end
        maximum(diam.(region)) < tol && break
    end
    return region
end

tfullcontract(regions, contractor, tol = 1e-7) =
    qmap(r -> fullcontract(r, contractor, tol), collect(regions))
