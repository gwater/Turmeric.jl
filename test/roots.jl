
using NumberIntervals, IntervalRootFinding2, StaticArrays, ForwardDiff
using Test

import IntervalRootFinding2: interval

function all_unique(rts)
    all(root_status.(rts) .== :unique)
end

function roots_dist(rt1::Root{T}, rt2::Root{T}) where {T}
    d = dist.(interval(rt1), interval(rt2))
    return sum(d)
end

function roots_dist(rt1::Root{T}, rt2::Root{T}) where {T <: Complex}
    dreal = dist(real(interval(rt1)), real(interval(rt2)))
    dimag = dist(imag(interval(rt1)), imag(interval(rt2)))

    return sum(dreal) + sum(dimag)
end

function test_newtonlike(f, deriv, X, method, nsol, tol=1e-10)
    rts = roots(f, X, method)
    @test length(rts) == nsol
    @test all_unique(rts)
    @test sum(roots_dist.(rts, roots(f, deriv, X, method))) < tol
end

newtonlike_methods = [Newton, Krawczyk]

@testset "1D roots" begin
    # Default
    rts = roots(sin, NumberInterval(-5, 5))
    @test length(rts) == 3
    @test all_unique(rts)

    # Bisection
    rts = roots(sin, NumberInterval(-5, 6), Bisection, 1e-3)
    @test length(rts) == 3

    # Refinement
    rts = roots(sin, rts, Newton)
    @test all_unique(rts)

    for method in newtonlike_methods
        test_newtonlike(sin, cos, NumberInterval(-5, 5), method, 3)
    end

    # Infinite interval
    rts = roots(x -> x^2 - 2, NumberInterval(-Inf, Inf))
    @test_broken length(rts) == 2
end


@testset "2D roots" begin
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)
    X = SVector(
        NumberInterval(-6, 6),
        NumberInterval(-6, 6)
    )

    # Bisection
    rts = roots(f, X, Bisection, 1e-3)
    @test length(rts) == 4

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(f, xx)
        test_newtonlike(f, deriv, X, method, 2)
    end

    # Infinite interval
    X = SVector(NumberInterval(-Inf, Inf), NumberInterval(-Inf, Inf))
    rts = roots(f, X, Newton)
    @test_broken length(rts) == 2
end


# From R docs:
# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/broyden

@testset "3D roots" begin
    function g(x)
        (x1, x2, x3) = x
        SVector(    x1^2 + x2^2 + x3^2 - 1,
                    x1^2 + x3^2 - 0.25,
                    x1^2 + x2^2 - 4x3
                )
    end

    X = NumberInterval(-5, 5)
    XX = SVector(X, X, X)

    for method in newtonlike_methods
        rts = roots(g, XX, method)
        @test length(rts) == 4
        @test all_unique(rts)
        deriv = xx -> ForwardDiff.jacobian(g, xx)
        @test rts == roots(g, deriv, XX, method)
    end
end

@testset "Out of domain" begin
    for method in newtonlike_methods
        @test length(roots(log, NumberInterval(-100, 2), method)) == 1
        @test length(roots(log, NumberInterval(-100, -1), method)) == 0
    end
end

@testset "Infinite domain" begin
    for method in newtonlike_methods
        rts = roots(x -> x^2 - 2, NumberInterval(-Inf, Inf), method)
        @test_broken length(filter(isunique, rts)) == 2
    end
end

@testset "NaN return value" begin
    f(xx) = ( (x, y) = xx; SVector(log(y/x) + 3x, y - 2x) )
    X = SVector(NumberInterval(-100, 100), NumberInterval(-100, 100))
    for method in newtonlike_methods
        rts = roots(f, X, method)
        @test length(filter(isunique, rts)) == 1
        @test length(filter(x -> contains_zero(x.interval), rts)) == 1
    end
end

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf = ∇(f)
    XX = SVector(NumberInterval(-5, 6), NumberInterval(-5, 6))
    tol = 1e-5

    for method in newtonlike_methods
        deriv = xx -> ForwardDiff.jacobian(gradf, xx)
        test_newtonlike(gradf, deriv, XX, method, 25, tol)
    end
end

@testset "Complex roots" begin
    X = NumberInterval(-5, 5)
    Xc = Complex(X, X)
    f(z) = z^3 - 1

    # Default
    rts = roots(f, Xc)
    @test length(rts) == 3

    # Bisection
    rts = roots(f, Xc, Bisection, 1e-3)
    @test length(rts) == 7

    for method in newtonlike_methods
        deriv = z -> 3*z^2
        test_newtonlike(f, deriv, Xc, method, 3, 1e-10)
    end
end
