#=
pytest/util.jl
Test julia version vs python version
Assumes that sigpy package is installed.
=#

using MRIPulses: dinf
using Test: @test, @testset


# Function to call the python implementation
function py_dinf(d1=nothing, d2=nothing)
    if isnothing(d1) && isnothing(d2)
        args = ""
    elseif isnothing(d2)
        args = "$d1"
    else
        args = "$d1, $d2"
    end

#   cmd = `python3 -c "import util; print(util.dinf($args))"` # local util.py
    cmd = `python3 -c "import sigpy.mri.rf.util as util; print(util.dinf($args))"`

    # Ensure current directory is in PYTHONPATH
    env = copy(ENV)
    env["PYTHONPATH"] = string(get(env, "PYTHONPATH", ""), ":", @__DIR__)

    output = read(setenv(cmd, env), String)
    return parse(Float64, strip(output))
end


@testset "util.jl vs util.py" begin
    # Test cases for d1 and d2
    d_vals = [0.1, 0.01, 0.001, 0.0001]

    for d1 in d_vals
        for d2 in d_vals
            jl_res = dinf(d1, d2)
            py_res = py_dinf(d1, d2)

            @test jl_res ≈ py_res
            if !(jl_res ≈ py_res)
                println("Mismatch at d1=$d1, d2=$d2: Julia=$jl_res, Python=$py_res")
            end
        end
    end

    # Test default arguments
    @test dinf() ≈ py_dinf()
end
