using MRIPulses: MRIPulses
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_ambiguities(MRIPulses) # if isolation needed
    Aqua.test_all(MRIPulses; ambiguities = false)
end
