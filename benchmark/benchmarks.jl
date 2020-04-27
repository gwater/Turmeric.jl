include("../examples/smiley_examples.jl")
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

using BenchmarkTools
using ForwardDiff
using IntervalArithmetic
using IntervalRootFinding2
using StaticArrays
import Random

const SUITE = BenchmarkGroup()

Random.seed!(0)  # Seed the RNG to get consistent results
tol = 1e-15

include("dietmar_ratz_functions.jl")

S = SUITE["Smiley"] = BenchmarkGroup()
for example in (SmileyExample22, SmileyExample52, SmileyExample54) #, SmileyExample55)
    s = S[example.title] = BenchmarkGroup()
    for method in (Newton, Krawczyk)
        s[string(method)] = @benchmarkable roots($(example.f), $(example.region), $method, $tol)
    end
end

S = SUITE["Rastigrin stationary points"] = BenchmarkGroup()

# Rastrigin function:
const A = 10

f(x, y) = 2A + x^2 - A*cos(2π*x) + y^2 - A*cos(2π*y)
f(X) = f(X...)

ForwardDiff.gradient(f, X::IntervalBox) = ForwardDiff.gradient(f, X.v)
∇f = X -> ForwardDiff.gradient(f, X)

L = 5.0
X = SVector(-L..(L+1), -L..(L+1))

for method in (Newton, Krawczyk)
    S[string(method)] = @benchmarkable roots($(∇(f)), $X, $method, 1e-5)
end

S = SUITE["Dietmar-Ratz"] = BenchmarkGroup()
X = Interval(0.75, 1.75)

for (k, dr) in enumerate(dr_functions)
    s = S["Dietmar-Ratz $k"] = BenchmarkGroup()

    if k != 8  # dr8 is excluded as it has too many roots
        for method in (Newton, Krawczyk)
            s[string(method)] = @benchmarkable roots($dr, $X, $method, $tol)
        end
    end

    s["Automatic differentiation"] = @benchmarkable ForwardDiff.derivative($dr, $X)
end
