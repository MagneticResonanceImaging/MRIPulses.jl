#=
slr/matrix.jl
matrix helper functions
=#


"""
    toeplitz(c, r=c)

Construct a Toeplitz matrix from column `c` and row `r`.
"""
function toeplitz(c::AbstractVector, r::AbstractVector=c)
    nc = length(c)
    nr = length(r)
    res = zeros(eltype(c), nc, nr)
    for j in 1:nr, i in 1:nc
        if i ≥ j
            res[i, j] = c[i - j + 1]
        else
            res[i, j] = r[j - i + 1]
        end
    end
    return res
end


"""
    hankel(c, r=zeros(eltype(c), length(c)))

Construct a Hankel matrix from column `c` and row `r`.
"""
function hankel(c::AbstractVector, r::AbstractVector=zeros(eltype(c), length(c)))
    nc = length(c)
    nr = length(r)
    res = zeros(eltype(c), nc, nr)
    for j in 1:nr, i in 1:nc
        idx = i + j - 1
        if idx ≤ nc
            res[i, j] = c[idx]
        else
            res[i, j] = r[idx - nc + 1]
        end
    end
    return res
end


"""
    hadamard(n)

Generate a Hadamard matrix of size `n` (must be a power of 2).
"""
function hadamard(n::Int; T::Type{<:Number} = Float64)
    if n == 1
        return ones(T, 1, 1)
    end
    ispow2(n) || throw("n = $n not a power of 2")
    Hn2 = hadamard(n >> 1)
    return [Hn2 Hn2; Hn2 -Hn2]
end
