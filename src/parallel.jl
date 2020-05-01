using Base.Threads

export troots, trefine

struct BisectionLimit end

struct ThreadBuffer
    root_regions
    indeterminate_regions
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

    # CONTRACTABLE?
    if !contraction_empty
        region = region .∩ contraction
    end

    # TOLERANCE LIMIT?
    if diam(region) < tol
        push!(buffer.indeterminate_regions, region)
        return nothing
    end

    # RECURSION LIMIT?
    if generation > maxgenerations || Sys.free_memory() / 2^30 < 1
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

function troots(region, contractor, f, tol, maxgenerations = 20)
    # setup data structures
    buffers = map(i -> ThreadBuffer(region), 1:nthreads())

    push!(buffers[1].indeterminate_regions, region)
    finished = false
    # this we can turn into an iterator
    while !finished
        finished = true
        vec_regions = qmap(buffers) do buffer
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
            task = troots!(buffers, _region, contractor, f, tol, 1, maxgenerations)
        end
        # wait for cascade to finish
        finished = isfinished(tasks)
        # repeats if we have unfinished regions
    end

    # collect results
    root_regions = reduce(vcat, map(b -> b.root_regions, buffers))
    indeterminate_regions =
        reduce(vcat, map(b -> b.indeterminate_regions, buffers))
    return root_regions, indeterminate_regions
end

function trefine(region, contractor, tol, maxiters = 100)
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
