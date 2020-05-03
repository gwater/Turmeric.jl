# This file is MIT licensed

module IntervalRootFinding2

using ForwardDiff
using StaticArrays

using LinearAlgebra: I, Diagonal

import Base: ⊆, show, big, \

## Root finding
export roots, bisect

import IntervalArithmetic: mid, ±, isdisjoint, isinterior, entireinterval,
    where_bisect, diam, bisect

diam(a::AbstractVector) = maximum(diam.(a))

function bisect(X::AbstractVector, α=where_bisect)
    i = argmax(diam.(X))  # find longest side
    return bisect(X, i, α)
end

function bisect(X::AbstractVector, i::Integer, α=where_bisect)
    x1, x2 = bisect(X[i], α)
    X1 = setindex(X, x1, i)
    X2 = setindex(X, x2, i)
    return (X1, X2)
end

contains_root(region) = all(0 .∈ region)
isdisjoint(a::AbstractVector, b::AbstractVector) = any(isdisjoint.(a, b))
isinterior(a::AbstractVector, b::AbstractVector) = all(isinterior.(a, b))

strict_isinterior(a, b) = isinterior(a, b) && !any(a .=== entireinterval.(a))

include("contractors.jl")

const default_tolerance = 1e-7
const default_contractor = Krawczyk()

include("parallel.jl")
include("roots.jl")

end
