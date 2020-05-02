# This file is part of the ValidatedNumerics.jl package; MIT licensed

module IntervalRootFinding2

using ForwardDiff
using StaticArrays

using LinearAlgebra: I, Diagonal

import Base: ⊆, show, big, \

## Root finding
export
    derivative, jacobian,  # reexport derivative from ForwardDiff
    Root, is_unique,
    roots, find_roots,
    bisect

export isunique, root_status

const derivative = ForwardDiff.derivative
const D = derivative

import IntervalArithmetic: where_bisect

include("root_object.jl")

include("complex.jl")
include("contractors.jl")
include("branch_and_bound.jl")
include("roots.jl")

include("parallel.jl")

gradient(f) = X -> ForwardDiff.gradient(f, X[:])

const ∇ = gradient

export ∇

function find_roots(f::Function, a::T, method::Function = newton;
                    tolerance = eps(T), debug = false, maxlevel = 30) where {T}

    method(f, a; tolerance=tolerance, debug=debug, maxlevel=maxlevel)
end

function find_roots(f::Function, f_prime::Function, a::T, method::Function=newton;
                    tolerance=eps(T), debug=false, maxlevel=30) where {T}

    method(f, f_prime, a; tolerance=tolerance, debug=debug, maxlevel=maxlevel)
end

function find_roots_midpoint(f::Function, region, method::Function=newton;
           tolerance=eps(1.0*a), debug=false, maxlevel=30, precision=-1)

    roots = find_roots(f, region, method; tolerance=tolerance, debug=debug, maxlevel=maxlevel, precision=precision)

    T = eltype(roots[1].interval)

    midpoints = T[]
    radii = T[]
    root_symbols = Symbol[]  # :unique or :unknown

    if length(roots) == 0
        return (midpoints, radii, root_symbols)  # still empty
    end

    for root in roots
        midpoint, radius = midpoint_radius(root.interval)
        push!(midpoints, midpoint)
        push!(radii, radius)

        push!(root_symbols, root.status)

    end

    (midpoints, radii, root_symbols)

end


function clean_roots(f, roots)

    # order, remove duplicates, and include intervals X only if f(X) contains 0
    sort!(roots, lt=lexless)
    roots = unique(roots)
    roots = filter(x -> 0 ∈ f(x.interval), roots)

    # merge neighbouring roots if they touch:

    if length(roots) < 2
        return roots
    end


    new_roots = eltype(roots)[]

    base_root = roots[1]

    for i in 2:length(roots)
        current_root = roots[i]

        if isempty(base_root.interval ∩ current_root.interval) ||
                (base_root.status != current_root.status)

            push!(new_roots, base_root)
            base_root = current_root
        else
            base_root = Root(hull(base_root.interval, current_root.interval), base_root.status)
        end
    end

    push!(new_roots, base_root)

    return new_roots

end

end
