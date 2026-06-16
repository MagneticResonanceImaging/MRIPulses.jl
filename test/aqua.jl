using MRIPulses: MRIPulses
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(MRIPulses)
end
