using LinearAlgebra


export Bisection, Newton, Krawczyk, contract

struct Newton end
struct Krawczyk end
struct Bisection end
const NewtonLike = Union{Type{Newton}, Type{Krawczyk}}

"""
    contract(f, f′, region, method, α)

Contract region `region` using `method` which can be Newton(), Krawczyk(),
or Bisection(). (Note: Bisection() is a no-op implemented for benchmarking.)
"""
contract(f, f′, region, ::Bisection, α) = region

function contract(f, f′, region::T, ::Newton, α) where {T}
    m = convert(T, mid(region, α))
    return m - (f(m) / f′(region))
end

function contract(f::Function, jacobian::Function, region::AbstractVector{T}, ::Newton, α) where T
    m = T.(mid.(region, α))
    J = jacobian(region)
    try
        return convert(typeof(region), m .- (J \ f(m)))
    catch
        return region .± Inf
    end
end


function contract(f, f′, region::T, ::Krawczyk, α) where T
    m = convert(T, mid(region, α))
    Y = 1 / f′(m)

    return m - Y*f(m) + (1 - Y*f′(region)) * (region - m)
end

function contract(f, jacobian, region::T, ::Krawczyk, α) where T <: AbstractVector
    m = mid.(region, α)
    mm = convert(T, m)
    J = jacobian(region)
    try
        Y = mid.(inv(jacobian(mm)))
        return m - Y*f(mm) + (I - Y*J) * (region - m)
    catch
        return region .± Inf
    end
end
