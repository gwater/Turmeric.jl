using LinearAlgebra

import IntervalArithmetic: mid, ±, isdisjoint

export Contractor
export Bisection, Newton, Krawczyk

"""
    𝒩(f, f′, X, α)

Single-variable Newton operator.

The symbol for the operator is accessed with `\\scrN<tab>`.
"""
function 𝒩(f, f′, X::T, α) where {T}
    m = convert(T, mid(X, α))
    return m - (f(m) / f′(X))
end

"""
    𝒩(f, jacobian, X, α)

Multi-variable Newton operator.
"""
function 𝒩(f::Function, jacobian::Function, X::AbstractVector{T}, α) where T # multidimensional Newton operator
    m = T.(mid.(X, α))
    J = jacobian(X)
    try
        return convert(typeof(X), m .- (J \ f(m)))
    catch
        return X .± Inf
    end
end

"""
    𝒦(f, f′, X, α)

Single-variable Krawczyk operator.

The symbol for the operator is accessed with `\\scrK<tab>`.
"""
function 𝒦(f, f′, X::T, α) where T
    m = convert(T, mid(X, α))
    Y = 1 / f′(m)

    return m - Y*f(m) + (1 - Y*f′(X)) * (X - m)
end

"""
    𝒦(f, jacobian, X, α)

Multi-variable Krawczyk operator.
"""
function 𝒦(f, jacobian, X::AbstractVector, α)
    m = mid.(X, α)
    mm = convert(typeof(X), m)
    J = jacobian(X)
    Y = mid.(inv(jacobian(mm)))

    return m - Y*f(mm) + (I - Y*J) * (X - m)
end


"""
    Contractor{F}

Abstract type for contractors.
"""
abstract type Contractor{F} end

"""
    Bisection{F} <: Contractor{F}

Contractor type for the bisection method.
"""
struct Bisection{F} <: Contractor{F}
    f::F
end

function (contractor::Bisection)(r, tol)
    X = interval(r)
    image = (contractor.f)(X)

    if root_status(r) == :empty || !all(0 .∈ image)
        return Root(X, :empty)
    end

    return Root(X, :unknown)
end

for (Method, 𝒪) in ((:Newton, 𝒩), (:Krawczyk, 𝒦))
    doc = """
        $Method{F, FP} <: Contractor{F}

    Contractor type for the $Method method.

    # Fields
      - `f::F`: function whose roots are searched
      - `f::FP`: derivative or jacobian of `f`

    -----

        (C::$Method)(X, tol, α=where_bisect)

    Contract an interval `X` using $Method operator and return the
    contracted interval together with its status.

    # Inputs
      - `R`: Root object containing the interval to contract.
      - `tol`: Precision to which unique solutions are refined.
      - `α`: Point of bisection of intervals.
    """

    @eval begin
        struct $Method{F, FP} <: Contractor{F}
            f::F
            f′::FP   # use \prime<TAB> for ′
        end

        function (C::$Method)(R, tol, α=where_bisect)
            op = x -> $𝒪(C.f, C.f′, x, α)
            rt = determine_region_status(op, C.f, R)
            return refine(op, rt, tol)
        end

        @doc $doc $Method
    end
end


"""
    safe_isempty(X)

Similar to `isempty` function for `IntervalBox`, but also works for `SVector`
of `Interval`.
"""
safe_isempty(X) = isempty(IntervalBox(X))

"""
    determine_region_status(contract, f, R)

Contraction operation for contractors using the first derivative of the
function.

Currently `Newton` and `Krawczyk` contractors use this.
"""
determine_region_status(op, f, R) =
    determine_region_status(op, f, interval(R), root_status(R))

const actions = (
    unknown = :bisect_contraction,
    image_nonzero = :discard,
    image_empty = :discard,
    empty_contraction = :discard,
    disjoint_contraction = :discard,
    interior_contraction = :mark_unique,
)

function assess_image(imX)
    if !all(0 .∈ imX)
        return :image_nonzero
    elseif any(isempty.(imX))
        return :image_empty
    end
    return :unknown
end

function assess_contraction(contracted_X, X)
    if any(isempty.(contracted_X))
        return :empty_contraction
    elseif any(isdisjoint.(contracted_X, X))
        return :disjoint_contraction
    elseif all(isinterior.(contracted_X, X))
        return :interior_contraction
    end
    return :unkown
end

# general rules:
# * always prefer discarding to bisecting
# * avoid op and f (potentially expensive)
# * refine only if f is defined everywhere on X (otherwise bisect)
# * it should never happen that we return empty interval as unique
# * no unique or empty state should be lost
# note:
# * how can f be undefined without exceptions? special definition?
# * why is this so extremely defensive? cant we exclude some cases?
function determine_region_status(op, f, X, former_status)
    if former_status in (:empty, :unique)
        # no further work required
        return Root(X, former_status)  # no-op
    end

    # refactor into assess_image()
    imX = f(X)
    if !all(0 .∈ imX)
        # no root inside X
        return Root(X, :empty)  # discard
    elseif any(isempty.(imX))
        # f is undefined on X
        return Root(X, :empty)  # discard
    end

    # not sure how I feel about this – might return empty set
    # also this happens at most once; we should just test this initially
    if any(isinf.(X))
        # force refine to constrain infinity
        return Root(op(X) .∩ X, :unknown)  # refined & bisect
    end

    contracted_X = op(X)

    # refactor into assess_contracted()
    if any(isempty.(contracted_X))
        # f is undefined on parts of X
        return Root(X, :unknown)  # bisect
    elseif any(isdisjoint.(contracted_X, X))
        # X and contracted_X are disjoint
        return Root(X, :empty)  # discard
    elseif all(isinterior.(contracted_X, X))
        # found unique root
        return Root(contracted_X, :unique)  # refined
    end

    # fallback – we've learned nothing
    return Root(contracted_X .∩ X, :unknown) # refined & bisect
end

"""
    refine(op, X, tol)

Generic refine operation for Krawczyk and Newton.
This function assumes that it is already known that `X` contains a unique root.
"""
function refine(op, X, tol)
    while maximum(diam.(X)) > tol  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = op(X) .∩ X
        all(NX .=== X) && break  # reached limit of precision
        any(isempty.(NX)) && break # dont exclude last element
        X = NX
    end

    return X
end

"""
    refine(op, X::Tuple{Symbol, Region}, tol)

Wrap the refine method to leave unchanged intervals that are not guaranteed to
contain an unique solution.
"""
function refine(op, R::Root, tol)
    root_status(R) != :unique && return R
    return Root(refine(op, interval(R), tol), :unique)
end
