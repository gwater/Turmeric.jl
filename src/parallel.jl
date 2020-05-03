using Base.Threads
using ThreadPools
using LazyArrays
using ForwardDiff

import Base: last, iterate
export last, iterate

export tfullcontract, fullcontract, ThreadedRootSearch

struct BisectionLimit end

struct ThreadBuffer{T, V <: AbstractVector{T}}
    root_regions::V
    indeterminate_regions::V
end
ThreadBuffer(region::T) where T = ThreadBuffer(T[], T[])

function troots!(buffers, region, _contract, f, tol, generation, maxgenerations)

    buffer = buffers[threadid()]

    # DISCARD?
    image = f(region)
    if !contains_root(image) || isempty(image)
        return nothing
    end
    contraction = _contract(region)
    contraction_empty = any(isempty.(contraction))
    if  !contraction_empty && isdisjoint(contraction, region)
        return nothing
    end

    # UNIQUE?
    if  !contraction_empty && strict_isinterior(contraction, region)
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
    if generation > maxgenerations
        push!(buffer.indeterminate_regions, region)
        return BisectionLimit()
    end

    # BISECT.
    a, b = bisect(region)
    task1 = Threads.@spawn troots!(buffers, a, _contract, f, tol, generation + 1, maxgenerations)
    task2 = Threads.@spawn troots!(buffers, b, _contract, f, tol, generation + 1, maxgenerations)
    return task1, task2
end

isfinished(::Nothing) = true
isfinished(::BisectionLimit) = false
isfinished(t::Task) = isfinished(fetch(t))
isfinished(ts::Tuple{Task, Task}) = isfinished(ts[1]) & isfinished(ts[2])
isfinished(ts::AbstractVector) = mapreduce(isfinished, &, ts)

struct ThreadedRootSearch{B}
    buffers::B
    _contract
    f
    tol::Float64
    target_task_number::Int
    function ThreadedRootSearch(
        f,
        region::T,
        _contract,
        tol = default_tolerance,
        target_task_number = nothing
    ) where T
        n = nthreads()
        buffers = map(i -> ThreadBuffer(region), Tuple(1:n))
        push!(buffers[1].indeterminate_regions, region)
        return new{typeof(buffers)}(
            buffers,
            _contract,
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
                search._contract,
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

function _troots(f, region, _contract, tol = default_tolerance, target_task_number = nothing)
    search = ThreadedRootSearch(f, region, _contract, tol, target_task_number)
    return last(search)
end

function troots(
    f,
    region,
    contractor = default_contractor,
    tol = default_tolerance,
    target_task_number = nothing
)
    function deriv(region)
        return ForwardDiff.derivative(f, region)
    end
    _contract(region) = contract(f, deriv, region, contractor, where_bisect)
    return _troots(f, region, _contract, tol, target_task_number)
end

function troots(
    f,
    region::T,
    contractor = default_contractor,
    tol = default_tolerance,
    target_task_number = nothing
) where T <: AbstractVector
    function grad(region)
        jac = similar(region, Size(length(region), length(region)))
        ForwardDiff.jacobian!(jac, f, region)
        return jac
    end
    _contract(region) = contract(f, grad, region, contractor, where_bisect)
    return _troots(f, region, _contract, tol, target_task_number)
end

function fullcontract(_contract, region::T, tol = default_tolerance, maxiters = 100) where T
    for i in 1:maxiters
        region2 = region .∩ _contract(region)
        if maximum(diam.(region2)) < maximum(diam.(region))
            region = region2
        else
            break
        end
        maximum(diam.(region)) < tol && break
    end
    return region
end

tfullcontract(_contract, regions, tol = default_tolerance) =
    qmap(r -> fullcontract(_contract, r, tol), collect(regions))
