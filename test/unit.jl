using Test
using NumberIntervals
using Turmeric
using Base.Threads
using StaticArrays

import Turmeric: strict_isinterior, contains_root

@testset "guarded_reduce()" begin
    nothing
end

@testset "explore_region!()" begin
    nothing
end

@testset "ThreadedRootSearch iterator interface" begin
    nothing
end

@testset "pop_unfinished!()" begin
    nothing
end

@testset "collect_regions()" begin
    nothing
end

@testset "estimate_generations()" begin
    nothing
end

@testset "chunks(regions, chunksize)" begin
    nothing
end

@testset "bisect()" begin
    x1 = NumberInterval(-1., 1.)
    x2 = NumberInterval(1., 10.)
    c = SVector(x1, x2)
    a, b = bisect(c)
    @test all(a .∪ b .=== c)
    @test !all(a .=== c)
end

@testset "contains_root()" begin
    x1 = NumberInterval(-1., 1.)
    x2 = NumberInterval(1., 10.)
    @test contains_root(x1)
    @test !contains_root(x2)
    @test contains_root(SVector(x1, x1))
    @test !contains_root(SVector(x1, x2))
end

@testset "isinterior()" begin
    x1 = NumberInterval(1. ,2.)
    x2 = NumberInterval(0., 5.)
    @test isinterior(SVector(x1, x1), SVector(x2, x2))
    @test !isinterior(SVector(x1, x1), SVector(x2, x1))
end

@testset "strict_isinterior()" begin
    x = entireinterval()
    @test !strict_isinterior(x, x)
end

@testset "diam()" begin
    x1 = NumberInterval(1. ,2.)
    x2 = NumberInterval(-3., 5.)
    @test diam(SVector(x1, x2)) >= max(diam(x1), diam(x2))
    x3 = NumberInterval(1)
    @test diam(SVector(x3, x3)) == 0
end

@testset "Krawczyk contractor" begin
    nothing
end

@testset "GradientContractor" begin
    x = NumberInterval(-10.0, -1.0)
    contractor = GradientContractor(log, Newton(), x)
    @test_throws DomainError contractor(x)
end

@testset "trivial contractor" begin
    nothing
end

@testset "fullcontract()" begin
    nothing
end
