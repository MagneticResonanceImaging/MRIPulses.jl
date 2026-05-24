#=
firls.jl
=#


"""
    firls(numtaps, bands, desired; weight=ones(length(desired)÷2))

FIR filter design using least-squares error minimization.
Translated from SciPy's `signal.firls`.

# Arguments
- `numtaps::Int`: Number of taps in the FIR filter. Must be odd.
- `bands::AbstractVector`: Monotonic non-decreasing sequence of band edges.
- `desired::AbstractVector`: Desired gain at the start and end of each band.
- `weight::AbstractVector`: Relative weighting for each band.

# Returns
- `Vector{Float64}`: Optimal FIR filter coefficients.
"""
function firls(numtaps::Int, bands::AbstractVector, desired::AbstractVector; weight=ones(nbands))
    (iseven(numtaps) || numtaps < 1) && error("numtaps must be odd and ≥ 1")
    M = (numtaps - 1) ÷ 2

    bands = bands[:]
    desired = desired[:]

    iseven(length(bands)) || error("bands must contain frequency pairs.")

    nbands = length(bands) ÷ 2

    # Set up the linear matrix equation to be solved, Qa = b
    # q(n) = 1/π ∫ W(ω) cos(nω) dω (over 0->π)
    n = 0:(numtaps-1)

    # Calculate q(n)
    q = zeros(numtaps)
    for i in 1:nbands
        f1, f2 = bands[2*i-1], bands[2*i]
        @. q += weight[i] * (f2 * sinc(n * f2) - f1 * sinc(n * f1))
    end

    Q1 = toeplitz(q[1:M+1])
    Q2 = hankel(q[1:M+1], q[M+1:end])
    Q = Q1 + Q2

    # b(n) = ∫ W(f) D(f) cos(nπf) df (over 0->1)
    b = zeros(M + 1)
    for i in 1:nbands
        f1, f2 = bands[2*i-1], bands[2*i]
        d1, d2 = desired[2*i-1], desired[2*i]

        m = (d2 - d1) / (f2 - f1)
        c = d1 - f1 * m

        # n=0: ∫ (mf+c) df = mf^2/2 + cf
        b[1] += weight[i] * (m * (f2^2 - f1^2) / 2 + c * (f2 - f1))

        for k in 1:M
            nk = k
            term = f2 * (m * f2 + c) * sinc(nk * f2) - f1 * (m * f1 + c) * sinc(nk * f1)
            term += m * (cos(nk * π * f2) - cos(nk * π * f1)) / (nk * π)^2
            b[k+1] += weight[i] * term
        end
    end

    # Solve Q*a = b
    a = Q \ b

    # SciPy reconstruction: h = [a_M, ..., a_1, 2*a_0, a_1, ..., a_M]
    # Note: SciPy's a_0 is Julia's a[1].
    coeffs = [reverse(a[2:end]); 2 * a[1]; a[2:end]]
    return coeffs
end
