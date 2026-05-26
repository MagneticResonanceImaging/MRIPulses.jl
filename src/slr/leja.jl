#=
slr/leja.jl
Polynomial root ordering.
=#


"""
    leja(r)

Order roots `r` using Leja ordering for numerical stability.

See Reichel, LAA, 1991:
https://doi.org/10.1016/0024-3795(91)90386-B
"""
function leja(r::AbstractVector{T}) where T
    n = length(r)
    if n == 0 return r end
    r_leja = zeros(T, n)

    # Start with the root with the largest magnitude
    idx = argmax(abs.(r))
    r_leja[1] = r[idx]

    # Keep track of remaining roots
    mask = trues(n)
    mask[idx] = false

    for i in 2:n
        best_val = -Inf
        best_idx = -1
        for j in 1:n
            if mask[j]
                # val = product of distances to already selected roots
                # use log-sum for stability
                val = 0.0
                for k in 1:i-1
                    val += log(abs(r_leja[k] - r[j]) + 1e-16)
                end
                if val > best_val
                    best_val = val
                    best_idx = j
                end
            end
        end
        r_leja[i] = r[best_idx]
        mask[best_idx] = false
    end
    return r_leja
end
