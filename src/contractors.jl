using LinearAlgebra

import IntervalArithmetic: mid, ±

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
function determine_region_status(op, f, R)
    X = interval(R)
    former_status = root_status(R)

    imX = f(X)

    if former_status == :empty || !all(0 .∈ imX)
        return Root(X, :empty)
    end

    any(isempty.(imX)) && return Root(X, :empty)  # X is fully outside of the domain of f

    contracted_X = op(X)

    # Only happens if X is partially out of the domain of f
    any(isempty.(contracted_X)) && return Root(X, :unknown)  # force bisection

    # given that have the Jacobian, can also do mean value form
    NX = contracted_X .∩ X

    any(isinf.(X)) && return Root(NX, :unknown)  # force bisection
    any(isempty.(NX)) && return Root(X, :empty)

    if former_status == :unique || all(isinterior.(NX, X))  # isinterior; know there's a unique root inside
        return Root(NX, :unique)
    end

    return Root(NX, :unknown)
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
