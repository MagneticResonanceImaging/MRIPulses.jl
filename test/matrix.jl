# test/matrix.jl

using MRIPulses: hankel, toeplitz, hadamard
using Test: @test, @testset, @test_throws, @inferred

@testset "Matrix helpers" begin
    @testset "toeplitz" begin
        c = [1, 2, 3]
        T = @inferred toeplitz(c)
        @test T == [1 2 3; 2 1 2; 3 2 1]

        r = [π, 2im, 3//4]
        @inferred toeplitz(c, r)
    end

    @testset "hankel" begin
        c = [1, 2, 3]
        H = @inferred hankel(c)
        @test H == [1 2 3; 2 3 0; 3 0 0]
    end

    @testset "hadamard" begin
        @test @inferred hadamard(1) == [1.0;;]
        @test @inferred hadamard(2) == [1.0 1.0; 1.0 -1.0]
        @test_throws String hadamard(3)
    end
end
