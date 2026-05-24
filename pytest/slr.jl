#=
pytest/slr.jl
Test julia version vs python version.
Assumes that sigpy package is installed.
=#

using LinearAlgebra: norm
using MRIPulses: dinf, dzlp, dzls, dzmp, dzrf, msinc
using MRIPulses: b2a, b2rf, mag2mp, ab2rf, root_flip, leja
using MRIPulses: dz_gslider_rf, dz_hadamard_b, dz_recursive_rf, calc_ripples
#using MRIPulses: hankel, toeplitz, hadamard_mat

using Test
using JSON


# Helper to serialize Julia objects for JSON transfer to Python
function to_py_json(obj)
    if obj isa Complex
        return Dict("__complex__" => true, "real" => real(obj), "imag" => imag(obj))
    elseif obj isa AbstractArray
        return [to_py_json(x) for x in obj]
    elseif obj isa Dict
        return Dict(k => to_py_json(v) for (k, v) in obj)
    else
        return obj
    end
end


# Helper to call Python functions and return the result as a Julia object
function py_call_func(func_name, args...; kwargs...)
    # Serialize arguments to JSON strings
    json_args = replace(JSON.json([to_py_json(a) for a in args]), "'" => "\\'")
    json_kwargs = replace(JSON.json(Dict(string(k) => to_py_json(v) for (k, v) in kwargs)), "'" => "\\'")

    script = """
import sys
import numpy as np
import json
import os
import sigpy.mri.rf.slr as slr

def from_json(data):
    if isinstance(data, dict) and data.get('__complex__'):
        return complex(data['real'], data['imag'])
    if isinstance(data, list):
        return [from_json(x) for x in data]
    if isinstance(data, dict):
        return {k: from_json(v) for k, v in data.items()}
    return data

def to_numpy(obj):
    if isinstance(obj, list):
        return np.array(obj)
    return obj

args = [to_numpy(from_json(x)) for x in json.loads('$json_args')]
kwargs = {k: to_numpy(from_json(v)) for k, v in json.loads('$json_kwargs').items()}

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if np.iscomplexobj(obj):
                return {'__complex__': True, 'real': obj.real.tolist(), 'imag': obj.imag.tolist()}
            return obj.tolist()
        if isinstance(obj, (complex, np.complex128)):
            return {'__complex__': True, 'real': float(obj.real), 'imag': float(obj.imag)}
        return json.JSONEncoder.default(self, obj)

res = slr.$(func_name)(*args, **kwargs)
if isinstance(res, tuple):
    res = list(res)
print(json.dumps(res, cls=NumpyEncoder))
"""
    cmd = `python3 -c "$script"`
    env = copy(ENV)
    output = read(setenv(cmd, env), String)
    return reconstruct_complex(JSON.parse(output, dicttype=Dict))
end


function reconstruct_complex(data)
    if data isa Dict && get(data, "__complex__", false)
        return reconstruct_complex_parts(data["real"], data["imag"])
    elseif data isa AbstractArray
        return [reconstruct_complex(x) for x in data]
    elseif data isa Dict
        return Dict(k => reconstruct_complex(v) for (k, v) in data)
    else
        return data
    end
end


function reconstruct_complex_parts(r, i)
    if r isa AbstractArray && i isa AbstractArray
        return reconstruct_complex_parts.(r, i)
    else
        return Complex(Float64(r), Float64(i))
    end
end


function py_to_jl(res)
    if res isa AbstractArray
        if isempty(res)
            return res
        end
        if all(x -> x isa Complex, res)
            return ComplexF64.(res)
        elseif all(x -> x isa Number, res)
            return Float64.(res)
        elseif all(x -> x isa AbstractArray, res)
            if all(x -> x isa AbstractArray && length(x) == length(res[1]), res)
                return hcat(py_to_jl.(res)...)'
            else
                return [py_to_jl(x) for x in res]
            end
        else
            return res
        end
    else
        return res
    end
end


py_call_complex(args...) =
    py_to_jl(reconstruct_complex(py_call_func(args...))) # helper


@testset "calc_ripples" begin
    for ptype in [:st, :ex, :se, :inv, :sat]
        jl_res = calc_ripples(ptype)
        py_res = py_call_func("calc_ripples", string(ptype))
        @test [jl_res...] ≈ [py_res...]
    end
end


#=
sigpy msinc() uses a lazy "+ 0.00001" kludge, whereas julia uses a precise
version of sinc(), leading to differences necessitating a larger tolerance here.
=#
@testset "msinc" begin
    n, m = 64, 2.0
    jl_res = msinc(n, m)
    py_res = py_call_func("msinc", n, m)
    @test jl_res ≈ py_res atol=5e-6 # see note above
end


@testset "dzlp" begin
    n, tb = 64, 4
    jl_res = dzlp(n, tb)
    py_res = py_call_func("dzlp", n, tb)
    @test jl_res ≈ py_res
end


# Slight differences in remez() require slightly larger tolerance here
@testset "dzmp" begin
    n, tb = 64, 4
    jl_res = dzmp(n, tb)
    py_res = py_call_func("dzmp", n, tb)
    @test jl_res ≈ py_res atol=5e-6
end


@testset "dzls" begin
    n, tb = 64, 4
    jl_res = dzls(n, tb)
    py_res = py_call_func("dzls", n, tb)
    @test jl_res ≈ py_res atol=5e-6
end


@testset "mag2mp" begin
    x = abs.(randn(64))
    jl_res = mag2mp(x)
    py_res = py_call_func("mag2mp", x)
    @test jl_res ≈ py_res atol=5e-6
end


@testset "b2a" begin
    b = msinc(64, 1.0)
    jl_res = b2a(b)
    py_res = py_call_func("b2a", b)
    @test jl_res ≈ py_res atol=5e-6
end


@testset "ab2rf" begin
    b = msinc(64, 1.0) .* 0.1
    a = b2a(b)
    jl_res = ab2rf(a, b)
    py_res = py_call_func("ab2rf", a, b)
    @test jl_res ≈ py_res
end


@testset "dzrf" begin
    n = 64
    tb = 4

    test_cases = [
        (:st, :pm),
        (:st, :ms),
        (:st, :ls),
        (:st, :min),
        (:st, :max),
        (:ex, :ms),
        (:se, :ls),
        (:inv, :pm),
        (:sat, :min)
    ]

    for (pt, ft) in test_cases
        @testset "$pt $ft" begin
            jl_res = dzrf(n=n, tb=tb, ptype=pt, ftype=ft)
            py_res = py_call_func("dzrf", n, tb, string(pt), string(ft))
            # @show norm(jl_res - py_res) / norm(py_res)
            @test jl_res ≈ py_res rtol=3e-5 # Relaxed to account for norm accumulation
        end
    end
end


@testset "root_flip" begin
    # Use dzls to get a pulse with candidate roots
    n = 32
    tb = 4.0
    d1 = 0.01
    d2 = 0.01
    flip = π/2
    b = dzls(n, tb, d1, d2)

    jl_res_rf, jl_res_b = root_flip(b, d1, flip, tb)
    py_res = py_call_func("root_flip", b, d1, flip, tb)
    py_res_rf = py_to_jl(reconstruct_complex(py_res[1]))
    py_res_b = py_res[2]

    @test jl_res_rf ≈ real(jl_res_rf)
    @test py_res_rf ≈ real(py_res_rf)
    jl_res_rf = real(jl_res_rf)
    py_res_rf = real(py_res_rf)

    # root_flip() can find equivalent solutions with different orders
    # so reverse() could be needed here
    @test isapprox(jl_res_rf, py_res_rf, atol=5e-6) ||
          isapprox(jl_res_rf, reverse(py_res_rf), atol=5e-6)
    @test isapprox(jl_res_b, py_res_b, atol=5e-6) ||
          isapprox(jl_res_b, reverse(py_res_b), atol=5e-6)
end


@testset "dz_recursive_rf" begin
    # Design a recursive RF pulse
    n_seg = 2 # Smaller number for more stable comparison
    tb = 4
    n = 32
    jl_res = dz_recursive_rf(n_seg=n_seg, tb=tb, n=n)
    py_res = py_call_complex("dz_recursive_rf", n_seg, tb, n)

    @test jl_res ≈ real(jl_res)
    @test py_res ≈ real(py_res)
    jl_res = real(jl_res)
    py_res = real(py_res)

    @test jl_res ≈ py_res atol=5e-6
end
