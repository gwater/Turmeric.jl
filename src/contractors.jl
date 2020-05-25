using LinearAlgebra
using LinearAlgebra: I

using ForwardDiff
using StaticArrays

import IntervalArithmetic: where_bisect

export Bisection, Newton, Krawczyk, GradientContractor, TrivialContractor

function _newton_contract(f, derivative, region::T, mid_point) where {T}
    m = convert(T, mid(region, mid_point))
    return m - (f(m) / derivative(region))
end

function _newton_contract(f, jacobian, region::R, mid_point) where {T, R <: AbstractVector{T}}
    m = T.(mid.(region, mid_point))
    J = jacobian(region)
    fm = f(m)
    try
        return convert(R, m .- (J \ fm))
    catch
        return region .± Inf
    end
end

function _krawczyk_contract(f, derivative, region::T, mid_point) where T
    m = convert(T, mid(region, mid_point))
    Y = 1 / derivative(m)
    return m - Y * f(m) + (1 - Y * derivative(region)) * (region - m)
end

function _krawczyk_contract(f, jacobian, region::T, mid_point) where T <: AbstractVector
    m = mid.(region, mid_point)
    mm = convert(T, m)
    J = jacobian(region)
    fmm = f(mm)
    try
        Y = mid.(inv(jacobian(mm)))
        return m - Y * fmm + (I - Y * J) * (region - m)
    catch
        return region .± Inf
    end
end

struct Newton end
struct Krawczyk end
struct Bisection end

"""
    GradientContractor(f, method, jacobian)
    GradientContractor(f, method, region)

Contractor for `f` based on `method` (`Newton()` or `Krawczyk()`), either using
the explicitely supplied function `jacobian(region)` or a jacobian method for
`region` using automatic differentiation of `f` from `ForwardDiff.jl`.

    (contractor::GradientContractor)(region, mid_point = IntervalArithmetic.where_bisect)

Calling instances of `GradientContractor` returns a contraction of `region`.
"""
struct GradientContractor{T}
    f
    method::T
    jacobian
end

function GradientContractor(f, method, ::R) where R <: AbstractVector
    function grad(region)
        jac = similar(region, Size(length(region), length(region)))
        ForwardDiff.jacobian!(jac, f, region)
        return jac
    end
    return GradientContractor(f, method, grad)
end

function GradientContractor(f, method, ::R) where R <: Number
    derivative(region) = ForwardDiff.derivative(f, region)
    return GradientContractor(f, method, derivative)
end

(contractor::GradientContractor{Newton})(region, mid_point = where_bisect) =
    _newton_contract(contractor.f, contractor.jacobian, region, mid_point)

(contractor::GradientContractor{Krawczyk})(region, mid_point = where_bisect) =
    _krawczyk_contract(contractor.f, contractor.jacobian, region, mid_point)

"""
    TrivialContractor()

Trivial contractor used for benchmarking which solely relies on bisections.

    (contractor::TrivialContractor)(region) = region

Calling instances of `TrivialContractor` returns `region` itself.
"""
struct TrivialContractor end

(contractor::TrivialContractor)(region) = region
