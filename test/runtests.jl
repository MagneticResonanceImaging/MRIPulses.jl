# runtests.jl

using MRIPulses: MRIPulses
using Test: @test, @testset, detect_ambiguities

include("aqua.jl")

# todo

@testset "ambiguities" begin
    @test isempty(detect_ambiguities(MRIPulses))
end
