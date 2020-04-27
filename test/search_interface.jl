using NumberIntervals
using IntervalRootFinding2

@testset "Root search iterator interface" begin
    contractor = Newton(sin, cos)
    search = BreadthFirstSearch(
        NumberInterval(-10, 10),
        contractor,
        1e-3
    )

    # First iteration of the search
    elem, state = iterate(search)

    # Optional iterator methods
    @test eltype(search) != Any
end


@testset "Search strategies" begin
    # DepthFristSearch and BreadthFristSearch should return the same results
    X = NumberInterval(-5, 5)
    rts = roots(sin, cos, X, Newton, DepthFirstSearch)
    @test Set(rts) == Set(roots(sin, cos, X, Newton, BreadthFirstSearch))
end
