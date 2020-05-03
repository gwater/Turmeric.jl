
using NumberIntervals, IntervalRootFinding2, StaticArrays, ForwardDiff
using Test

all_unique(rts) = length(rts[2]) == 0

function test_newtonlike(f, X, method, nsol)
    rts = roots(f, X, method)
    @test length(rts[1]) == nsol
    @test all_unique(rts)
end

newtonlike_methods = (Newton(), Krawczyk())

@testset "1D roots" begin
    # Default
    rts = roots(sin, NumberInterval(-5, 5))
    @test length(rts[1]) == 3
    @test all_unique(rts)

    # Bisection()
    rts = roots(sin, NumberInterval(-5, 6), Bisection(), 1e-3)
    @test length(rts[2]) == 3

    for method in newtonlike_methods
        test_newtonlike(sin, NumberInterval(-5, 5), method, 3)
    end

    # Infinite interval
    rts = roots(x -> x^2 - 2, NumberInterval(-Inf, Inf), Newton())
    @test length(rts[1]) == 2
end

@testset "2D roots" begin
    f(x, y) = SVector(x^2 + y^2 - 1, y - 2x)
    f(X) = f(X...)
    X = SVector(
        NumberInterval(-6, 6),
        NumberInterval(-6, 6)
    )

    # Bisection()
    rts = roots(f, X, Bisection(), 1e-3)
    @test_broken length(rts[2]) == 2

    for method in newtonlike_methods
        test_newtonlike(f, X, method, 2)
    end

    # Infinite interval
    X = SVector(NumberInterval(-Inf, Inf), NumberInterval(-Inf, Inf))
    rts = roots(f, X, Newton())
    @test length(rts[1]) == 2
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
        @test length(rts[1]) == 4
        @test all_unique(rts)
    end
end

@testset "Out of domain" begin
    for method in newtonlike_methods
        @test length(roots(log, NumberInterval(-100, 2), method)[1]) == 1
        @test length(roots(log, NumberInterval(-100, -1), method)[1]) == 0
    end
end

@testset "Infinite domain" begin
    for method in newtonlike_methods
        rts = roots(x -> x^2 - 2, NumberInterval(-Inf, Inf), method)
        @test length(rts[1]) == 2
        @test all_unique(rts)
    end
end

@testset "NaN return value" begin
    f(xx) = ( (x, y) = xx; SVector(log(y/x) + 3x, y - 2x) )
    X = SVector(NumberInterval(-100, 100), NumberInterval(-100, 100))
    for method in newtonlike_methods
        rts = roots(f, X, method)
        @test length(rts[1]) == 1
        @test length(rts[2]) > 0 # not broken, problem is ill-defined
    end
end

@testset "Stationary points" begin
    f(xx) = ( (x, y) = xx; sin(x) * sin(y) )
    gradf(x) = ForwardDiff.gradient(f, x)
    XX = SVector(NumberInterval(-5, 6), NumberInterval(-5, 6))
    tol = 1e-5

    for method in newtonlike_methods
        test_newtonlike(gradf, XX, method, 25)
    end
end
