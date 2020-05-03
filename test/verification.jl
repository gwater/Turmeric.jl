include("../examples/smiley_examples.jl")

using Test

using NumberIntervals, IntervalRootFinding2
using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

const tol = 1e-6

#NOTE: Bisection method performs badly in all examples
for method in (Krawczyk(), Newton())

@testset "$method: $(SmileyExample22.title)" begin
    region = NumberInterval.(SmileyExample22.region)
    roots_found, rest_regions = roots(SmileyExample22.f, region, method, tol)
    @test length(roots_found) == 8
    @test length(rest_regions) == 0
    # no reference data for roots given
end

for example in ifelse(method == Newton(),
    (SmileyExample52, SmileyExample54, ), # SmileyExample55 takes too long with Newton
    (SmileyExample52, SmileyExample54, SmileyExample55),
)
    @testset "$method: $(example.title)" begin
        region = NumberInterval.(example.region)
        roots_found, rest_regions = roots(example.f, region, method, tol)
        @test length(roots_found) == length(example.known_roots)
        @test length(rest_regions) == 0
        for rf in roots_found
            # check there is exactly one known root for each found root
            @test sum(!any(isempty.(rk .∩ rf)) for rk in example.known_roots) == 1
        end
    end
end

end # method loop
