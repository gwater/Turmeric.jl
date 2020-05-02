using IntervalRootFinding2
using Test
using Base.Threads

if nthreads() == 1
    include("branch_and_bound.jl")
    include("search_interface.jl")
    include("roots.jl")
    include("test_smiley.jl")
elseif nthreads() > 1
    include("parallel.jl")
end
