
import IntervalArithmetic: mid, ±, isdisjoint, isinterior, entireinterval,
    where_bisect, diam, bisect

export bisect

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
isinterior(a::AbstractVector, b::AbstractVector) = all(isinterior.(a, b))

strict_isinterior(a, b) = isinterior(a, b) && !any(a .=== entireinterval.(a))

# these definitions break the concept of Vectors as sets of their elements
# (rather than Cartesian boxes)
_isempty(a) = any(isempty.(a))
_isdisjoint(a, b) = any(isdisjoint.(a, b))
