# test/slr.jl

using FFTW: fft
using MRIPulses: dinf, dzlp, dzls, dzmp, dzrf, msinc
using MRIPulses: b2a, b2rf, mag2mp, ab2rf, root_flip, leja
using MRIPulses: dz_gslider_rf, dz_hadamard_b, dz_recursive_rf, calc_ripples
using MRIPulses: hankel, toeplitz, hadamard
using Test: @test, @testset, @test_throws, @inferred


@testset "SLRUtil" begin
    @testset "dinf" begin
        # Test default values
        @test @inferred(dinf()) ≈ 1.944048
        # Test specific values
        @test @inferred(dinf(0.01, 0.01)) ≈ 1.944048
        @test @inferred(dinf(0.001, 0.001)) > 1.944048
    end
end


@testset "Internal Helpers" begin
    @testset "toeplitz" begin
        c = [1, 2, 3]
        T = @inferred toeplitz(c)
        @test T == [1 2 3; 2 1 2; 3 2 1]
    end

    @testset "hankel" begin
        c = [1, 2, 3]
        H = @inferred hankel(c)
        @test H == [1 2 3; 2 3 0; 3 0 0]
    end

    @testset "leja" begin
        r = [1.0, 0.5, 0.1]
        rl = @inferred leja(r)
        @test length(rl) == 3
        @test sort(abs.(rl)) == sort(abs.(r))
    end

    @testset "hadamard" begin
        @test @inferred hadamard(1) == [1.0;;]
        @test @inferred hadamard(2) == [1.0 1.0; 1.0 -1.0]
        @test_throws String hadamard(3)
    end
end


@testset "SLR Core" begin
    n = 64
    tb = 4

    @testset "msinc" begin
        ms = @inferred msinc(n, tb/4)
        @test length(ms) == n
        @test eltype(ms) <: Real
    end

    @testset "dzlp" begin
        h = @inferred dzlp(n, tb, 0.01, 0.01)
        @test length(h) == n
        @test any(h .!= 0)
    end

    @testset "dzls" begin
        h = @inferred dzls(n, tb, 0.01, 0.01)
        @test length(h) == n
    end

    @testset "dzmp" begin
        h = @inferred dzmp(n, tb, 0.01, 0.01)
        @test length(h) == n
    end

    @testset "dzrf" begin
        # Small tip
        rf_st = dzrf(n=n, tb=tb, ptype=:st, ftype=:ms) # @NOTinferred
        @test length(rf_st) == n

        # Excitation
        rf_ex = dzrf(n=n, tb=tb, ptype=:ex, ftype=:ls) # @NOTinferred
        @test length(rf_ex) == n
        @test eltype(rf_ex) <: Complex

        # Test different filter types
        @test length(dzrf(n=n, tb=tb, ftype=:pm)) == n
        @test length(dzrf(n=n, tb=tb, ftype=:min)) == n
        @test length(dzrf(n=n, tb=tb, ftype=:max)) == n

        # Test spin-echo pulse type (ptype other than :st or :ex)
        rf_se_slr = dzrf(n=n, tb=tb, ptype=:se, ftype=:ms) # @NOTinferred
        @test length(rf_se_slr) == n

        # Test error for unknown filter
        @test_throws ErrorException dzrf(ftype=:unknown)
    end

    @testset "Polynomial Operations" begin
        b = dzrf(n=n, tb=tb, ptype=:st, ftype=:ms) # @NOTinferred

        @testset "b2a" begin
            a = @inferred b2a(b)
            @test length(a) == n
        end

        @testset "b2rf" begin
            rf = @inferred b2rf(b)
            @test length(rf) == n

            rf_cancel = @inferred b2rf(b, true)
            @test length(rf_cancel) == n
        end

        @testset "mag2mp" begin
            x = abs.(fft(b))
            mp = @inferred mag2mp(x)
            @test length(mp) == length(x)
        end
    end

    @testset "ab2rf" begin
        a = ones(ComplexF64, 10)
        b = zeros(ComplexF64, 10)
        b[5] = 0.1
        rf = @inferred ab2rf(a, b)
        @test length(rf) == 10
    end
end


@testset "Advanced Designs" begin
    @testset "root_flip" begin
        n = 32
        tb = 4
        b = @inferred msinc(n, tb/4)
        rf_out, b_out = root_flip(b, 0.01, π/2, tb, verbose=true) # verbose=true for coverage # @NOTinferred
        @test length(rf_out) >= n-2 # Allow for small zero-taps truncation
        @test length(b_out) >= n-2
    end

    @testset "dz_gslider_rf" begin
        # Test g=1
        rf1 = @inferred dz_gslider_rf(n=64, g=1, flip=π/2)
        @test size(rf1) == (64, 1)

        # Test g=5 (exercises non-centered logic in dz_gslider_b)
        rf5 = @inferred dz_gslider_rf(n=64, g=5, flip=π/2, tb=4)
        @test size(rf5) == (64, 5)
        @test any(rf5 .!= 0)
    end

    @testset "dz_hadamard_b" begin
        n = 64
        g = 4
        # gind=1
        b1 = dz_hadamard_b(n, g, 1, 4, 0.01, 0.01, 16) # @NOTinferred
        @test length(b1) == n
        # gind > 1
        b2 = dz_hadamard_b(n, g, 2, 4, 0.01, 0.01, 16)
        @test length(b2) == n
        @test any(b2 .!= 0)
    end
end


@testset "dz_recursive_rf" begin
    n_seg = 3
    tb = 4
    n = 32

    # Test gradient echo version, use_mz=true
    rf = dz_recursive_rf(n_seg=n_seg, tb=tb, n=n, z_pad_fact=2.0, use_mz=true)
    @test size(rf, 2) == n_seg

    # Test gradient echo version, use_mz=false
    rf_ideal = dz_recursive_rf(n_seg=n_seg, tb=tb, n=n, z_pad_fact=2.0, use_mz=false)
    @test size(rf_ideal, 2) == n_seg

    # Test spin echo version
    rf_se, rf_ref = dz_recursive_rf(n_seg=n_seg, tb=tb, n=n, se_seq=true, z_pad_fact=2.0)
    @test size(rf_se, 2) == n_seg
    @test length(rf_ref) > 0
end


@testset "Error Cases and Placeholders" begin
    @test_throws ErrorException calc_ripples(:unknown)
    @test_throws ErrorException dzrf(n=64, tb=4, ftype=:unknown)
    # Test all recognized pulse types for coverage
    @test length(calc_ripples(:st)) == 3
    @test length(calc_ripples(:ex)) == 3
    @test length(calc_ripples(:se)) == 3
    @test length(calc_ripples(:inv)) == 3
    @test length(calc_ripples(:sat)) == 3
end
