#=
matlab.jl
Test Julia versus Matlab version

Steps for comparing Matlab to Julia versions

1. git Toppe repo
git clone git@github.com:toppeMRI/toppe.git

2.
Launch Matlab

3.
Add top toppe directory to Matlab path, something like:
>> addpath ~/dat/git/clone/toppe

4.
Design an RF pulse with Matlab version and write to a file
>> slr6 = single(toppe.utils.rf.jpauly.dzrf(100, 6));
>> plot(slr6)
>> save slr6.mat slr6

5.
Run the test below.

For this test case, the maximum of mat6 is 0.059
and the maximum difference between the Julia and Matlab versions is 9.29e-5
which is a bit larger than I might have hoped
but if you plot them they are visually indistinguishable.
=#

using MAT: matread
using MRIPulses: dzrf
using Test: @test, @testset

@testset "matlab" begin
    mat6 = matread("slr6.mat")["slr6"][:,1]
    n = 100
    tb = 6
    julia6 = dzrf( ; n, tb, ptype = :st, ftype = :ls, )
    norm = x -> maximum(abs, x) # ∞ norm
    @test isapprox(julia6, mat6; atol = 1e-4, norm)
    # plot([julia6 mat6]) # visually indistinguishable
end
