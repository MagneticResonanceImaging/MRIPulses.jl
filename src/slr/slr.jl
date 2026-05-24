#=
slr.jl
Shinnar-Le Roux (SLR) RF pulse design for MRI
Initial translation from SigPy performed by Gemini 0.39.0.
https://github.com/mikgroup/sigpy/blob/main/sigpy/mri/rf/slr.py
=#

using DSP: remez, freqresp, PolynomialRatio, hilbert
using DSP.Windows: blackman
using FFTW: fft, ifft, ifftshift, fftshift
using Polynomials: Polynomial, roots, fromroots, coeffs

export dzrf, dzls, msinc, dzmp, fmp, dzlp, b2rf, b2a, mag2mp, ab2rf,
       dz_gslider_b, dz_gslider_rf, root_flip, dz_recursive_rf,
       dz_hadamard_b, calc_ripples


# --- Internal Helper Functions ---

"""
    sp_fft(x, center=true)

SigPy-style FFT. Defaults to centered FFT.
"""
function sp_fft(x::AbstractArray; center::Bool=true)
    if center
        return fftshift(fft(ifftshift(x)))
    else
        return fft(x)
    end
end

"""
    sp_ifft(x, center=true)

SigPy-style IFFT. Defaults to centered IFFT.
"""
function sp_ifft(x::AbstractArray; center::Bool=true)
    if center
        return fftshift(ifft(ifftshift(x)))
    else
        return ifft(x)
    end
end


# --- Primary Design Functions ---

"""
    dzrf(n=64, tb=4, ptype=:st, ftype=:ls, d1=0.01, d2=0.01, cancel_alpha_phs=false)

Primary function for design of pulses using the Shinnar-Le Roux (SLR) algorithm.

# Arguments
- `n::Int`: Number of time points.
- `tb::Real`: pulse time bandwidth product.
- `ptype::Symbol`: pulse type, `:st` (small-tip), `:ex` (π/2), `:se` (spin-echo), `:inv` (inversion), `:sat` (saturation).
- `ftype::Symbol`: type of filter: `:ms` (sinc), `:pm` (Parks-McClellan), `:min` (minphase), `:max` (maxphase), `:ls` (least squares).
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.
- `cancel_alpha_phs::Bool`: For `:ex` pulses, cancel alpha phase for a flatter profile.

# Returns
- `Vector{ComplexF64}`: Designed RF pulse.
"""
function dzrf(;
    n::Int = 64,
    tb::Real = 4,
    ptype::Symbol = :st,
    ftype::Symbol = :ls,
    d1::Real = 0.01,
    d2::Real = 0.01,
    cancel_alpha_phs::Bool = false,
)
    bsf, d1, d2 = calc_ripples(ptype, d1, d2)

    if ftype === :ms
        b = msinc(n, tb / 4)
    elseif ftype === :pm
        b = dzlp(n, tb, d1, d2)
    elseif ftype === :min
        b = reverse(dzmp(n, tb, d1, d2))
    elseif ftype === :max
        b = dzmp(n, tb, d1, d2)
    elseif ftype === :ls
        b = dzls(n, tb, d1, d2)
    else
        error("Filter type ($ftype) is not recognized.")
    end

    if ptype === :st
        rf = b
    elseif ptype === :ex
        rf = b2rf(bsf .* b, cancel_alpha_phs)
    else
        rf = b2rf(bsf .* b)
    end

    return rf
end


"""
    calc_ripples(ptype=:st, d1=0.01, d2=0.01)

Calculate effective SLR ripple levels for a specific pulse type.

# Arguments
- `ptype::Symbol`: Pulse type.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.

# Returns
- `Tuple{Real, Real, Real}`: `(bsf, d1, d2)`.
"""
function calc_ripples(ptype::Symbol = :st, d1::Real = 0.01, d2::Real = 0.01)
    if ptype === :st
        bsf = 1.0
    elseif ptype === :ex
        bsf = sqrt(0.5)
        d1 = sqrt(d1 / 2)
        d2 = d2 / sqrt(2)
    elseif ptype === :se
        bsf = 1.0
        d1 /= 4
        d2 = sqrt(d2)
    elseif ptype === :inv
        bsf = 1.0
        d1 /= 8
        d2 = sqrt(d2 / 2)
    elseif ptype === :sat
        bsf = sqrt(0.5)
        d1 /= 2
        d2 = sqrt(d2)
    else
        error("Pulse type ($ptype) is not recognized.")
    end
    return bsf, d1, d2
end


# --- Support Functions ---

"""
    dzls(n=64, tb=4, d1=0.01, d2=0.01)

Design a linear-phase FIR filter using the least-squares error minimization method.

# Arguments
- `n::Int`: Number of time points.
- `tb::Real`: Time-bandwidth product.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.

# Returns
- `Vector{Float64}`: Designed filter coefficients (truncated to `n` samples).
"""
function dzls(n::Int = 64, tb::Real = 4, d1::Real = 0.01, d2::Real = 0.01)
    di = dinf(d1, d2)
    w = di / tb
    f_edges = [0.0, (1 - w) * (tb / 2), (1 + w) * (tb / 2), (n / 2)] ./ (n / 2)

    h = firls(n + 1, f_edges, [1.0, 1.0, 0.0, 0.0]; weight=[1.0, d1 / d2])

    c = cis.(2π / (2 * (n + 1)) .* [0:n÷2; -n÷2:-1])
    h = real.(ifft(fft(h) .* c))

    # Lop off extra sample
    return h[1:n]
end


"""
    dzmp(n=64, tb=4, d1=0.01, d2=0.01)

Design a minimum-phase FIR filter.

# Arguments
- `n::Int`: Number of time points.
- `tb::Real`: Time-bandwidth product.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.

# Returns
- `Vector{ComplexF64}`: Designed minimum-phase filter coefficients.
"""
function dzmp(n::Int = 64, tb::Real = 4, d1::Real = 0.01, d2::Real = 0.01)
    n2 = 2 * n - 1
    di = 0.5 * dinf(2 * d1, 0.5 * d2^2)
    w = di / tb
    f_edges = [0.0, (1 - w) * (tb / 2), (1 + w) * (tb / 2), (n / 2)] ./ n
    weight = [1.0, 2 * d1 / (0.5 * d2^2)]
    hl = remez(n2, f_edges, [1.0, 0.0]; weight)
    return fmp(hl)
end


"""
    fmp(h)

Convert a linear-phase filter `h` to its minimum-phase equivalent.

# Arguments
- `h::AbstractVector`: Linear-phase filter coefficients.

# Returns
- `Vector{ComplexF64}`: Minimum-phase filter coefficients.
"""
function fmp(h::AbstractVector)
    ll = length(h)
    lp = Int(128 * 2^ceil(log2(ll)))
    pad_total = lp - ll
    hp = [zeros(Int(ceil(pad_total / 2))); h; zeros(Int(floor(pad_total / 2)))]
    # Match SigPy sp.fft(hp, center=True, norm=None)
    hpf = sp_fft(hp)
    hpfs = hpf .- minimum(real.(hpf)) * 1.000001
    hpfmp = mag2mp(sqrt.(abs.(hpfs)))
    # Match SigPy sp.ifft(ifftshift(conj(hpfmp)), center=False, norm=None)
    hpmp = ifft(ifftshift(conj.(hpfmp)))
    return hpmp[1:(ll + 1) ÷ 2]
end


"""
    dzlp(n=64, tb=4, d1=0.01, d2=0.01)

Design a linear-phase FIR filter using the Remez exchange algorithm.

# Arguments
- `n::Int`: Number of time points.
- `tb::Real`: Time-bandwidth product.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.

# Returns
- `Vector{Float64}`: Designed filter coefficients.
"""
function dzlp(n::Int = 64, tb::Real = 4, d1::Real = 0.01, d2::Real = 0.01)
    di = dinf(d1, d2)
    w = di / tb
    f_edges = [0.0, (1 - w) * (tb / 2), (1 + w) * (tb / 2), (n / 2)] ./ n
    weight = [1.0, d1 / d2]
    return remez(n, f_edges, [1.0, 0.0]; weight)
end


"""
    msinc(n=64, m=1.0)

Generate a Hamming-windowed sinc pulse with `m` lobes.

# Arguments
- `n::Int`: Number of time points.
- `m::Real`: Number of lobes.

# Returns
- `Vector{Float64}`: Windowed sinc pulse.
"""
function msinc(n::Int = 64, m::Real = 1.0)
    x = ((-n / 2):(n / 2 - 1)) / (n / 2)
    snc = @. sinc(2m * x)
    ms = @. snc * (0.54 + 0.46 * cos(π * x))
    return @. ms * 4 * m / n
end


"""
    dz_gslider_b(n=128, g=5, gind=1, tb=4, d1=0.01, d2=0.01, phi=π, shift=32)

Design a gSlider SLR beta parameter.

# Arguments
- `n::Int`: Number of time points.
- `g::Int`: Number of sub-slices.
- `gind::Int`: Sub-slice index.
- `tb::Real`: Time-bandwidth product.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.
- `phi::Real`: Sub-slice phase.
- `shift::Int`: Number of time points shift of pulse.

# Returns
- `Vector{ComplexF64}`: SLR beta parameter.

# References
- Setsompop, K. et al. (2018). 'High-resolution in vivo diffusion imaging of the
  human brain with generalized slice dithered enhanced resolution:
  Simultaneous multislice (gSlider-SMS)'. Magn. Reson. Med. 79, 141–151.
"""
function dz_gslider_b(n::Int = 128, g::Int = 5, gind::Int = 1,
     tb::Real = 4, d1::Real = 0.01, d2::Real = 0.01, phi::Real = π, shift::Int = 32,
)
    ftw = dinf(d1, d2) / tb
    c = cis.(2π / (2 * (n + 1)) .* [0:n÷2; -n÷2:-1])

    if isodd(g) && gind == (g + 1) ÷ 2
        if g == 1
            return dzls(n, tb, d1, d2)
        else
            f = [0, (1/g - ftw)*(tb/2), (1/g + ftw)*(tb/2), (1 - ftw)*(tb/2), (1 + ftw)*(tb/2), n/2] ./ (n/2)
            b_notch = firls(n + 1, f, [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]; weight=[1.0, 1.0, d1 / d2])
            b_sub = firls(n + 1, f, [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]; weight=[1.0, 1.0, d1 / d2])
            b = ifft(fft(b_notch .+ cis(phi) .* b_sub) .* c)
            return b[1:n]
        end
    else
        gcent = shift + (gind - g / 2 - 0.5) * tb / g
        if gind > 1 && gind < g
            f = [0, shift - (1 + ftw) * (tb / 2), shift - (1 - ftw) * (tb / 2),
                 gcent - (tb / g / 2 + ftw * (tb / 2)), gcent - (tb / g / 2 - ftw * (tb / 2)),
                 gcent + (tb / g / 2 - ftw * (tb / 2)), gcent + (tb / g / 2 + ftw * (tb / 2)),
                 shift + (1 - ftw) * (tb / 2), shift + (1 + ftw) * (tb / 2), (n / 2)]
            m_notch = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
            m_sub = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
            w = [d1 / d2, 1.0, 1.0, 1.0, d1 / d2]
        elseif gind == 1
            f = [0, shift - (1 + ftw) * (tb / 2), shift - (1 - ftw) * (tb / 2),
                 gcent + (tb / g / 2 - ftw * (tb / 2)), gcent + (tb / g / 2 + ftw * (tb / 2)),
                 shift + (1 - ftw) * (tb / 2), shift + (1 + ftw) * (tb / 2), (n / 2)]
            m_notch = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
            m_sub = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
            w = [d1 / d2, 1.0, 1.0, d1 / d2]
        else # gind == g
            f = [0, shift - (1 + ftw) * (tb / 2), shift - (1 - ftw) * (tb / 2),
                 gcent - (tb / g / 2 + ftw * (tb / 2)), gcent - (tb / g / 2 - ftw * (tb / 2)),
                 shift + (1 - ftw) * (tb / 2), shift + (1 + ftw) * (tb / 2), (n / 2)]
            m_notch = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
            m_sub = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
            w = [d1 / d2, 1.0, 1.0, d1 / d2]
        end
        f ./= (n / 2)

        b_n = hilbert(real.(ifft(fft(firls(n + 1, f, m_notch; weight=w)) .* c)[1:n]))
        b_s = hilbert(real.(ifft(fft(firls(n + 1, f, m_sub; weight=w)) .* c)[1:n]))
        b = b_n .+ cis(phi) .* b_s
        c_shift = @. cis(-2 * π / n * shift * (0:n-1)) / 2 * cis(-π / n * shift)
        return b .* c_shift
    end
end


"""
    dz_gslider_rf(n=256, g=5, flip=π/2, phi=π, tb=12, d1=0.01, d2=0.01, cancel_alpha_phs=true)

Design a gSlider RF pulse.
"""
function dz_gslider_rf(; n=256, g=5, flip=π/2, phi=π, tb=12, d1=0.01, d2=0.01, cancel_alpha_phs=true)
    bsf = sin(flip / 2)
    rf = zeros(ComplexF64, n, g)
    for gind in 1:g
        b = bsf .* dz_gslider_b(n, g, gind, tb, d1, d2, phi)
        rf[:, gind] = b2rf(b, cancel_alpha_phs)
    end
    return rf
end


"""
    dz_hadamard_b(n=128, g=5, gind=1, tb=4, d1=0.01, d2=0.01, shift=32)

Design a Hadamard-encoded SLR beta parameter.
"""
function dz_hadamard_b(n=128, g=5, gind=1, tb=4, d1=0.01, d2=0.01, shift=32)
    H = hadamard(g)
    encode = H[gind, :]
    ftw = dinf(d1, d2) / tb

    if gind == 1
        return dzls(n, tb, d1, d2)
    else
        f = [0.0, shift - (1 + ftw) * (tb / 2)]
        m = [0.0, 0.0]
        w = [d1 / d2]

        for ii in 1:g
            gcent = shift + (ii - g / 2 - 0.5) * tb / g
            if ii == 1 || encode[ii] != encode[ii-1]
                push!(f, gcent - (tb / g / 2 - ftw * (tb / 2)))
                push!(m, encode[ii])
            end
            if ii == g || encode[ii] != encode[ii+1]
                push!(f, gcent + (tb / g / 2 - ftw * (tb / 2)))
                push!(m, encode[ii])
                push!(w, 1.0)
            end
        end
        append!(f, [shift + (1 + ftw) * (tb / 2), n / 2])
        append!(m, [0.0, 0.0])
        push!(w, d1 / d2)
        f ./= (n / 2)

        c = cis.(2 * π / (2 * (n + 1)) .* [0:n÷2; -n÷2:-1])
        bp = firls(n + 1, f, Float64.(m .> 0); weight=w)
        bn = firls(n + 1, f, Float64.(m .< 0); weight=w)
        b = hilbert(real.(ifft(fft(bp .- bn) .* c)[1:n]))
        c_shift = @. cis(-2 * π / n * shift * (0:n-1)) / 2 * cis(-π / n * shift)
        return b .* c_shift
    end
end


"""
    b2rf(b, cancel_alpha_phs=false)

Convert an SLR beta parameter to an RF pulse.

# Arguments
- `b::AbstractVector`: SLR beta parameter.
- `cancel_alpha_phs::Bool`: Cancel alpha phase.

# Returns
- `Vector{ComplexF64}`: Designed RF pulse.
"""
function b2rf(b::AbstractVector, cancel_alpha_phs::Bool = false)
    a = b2a(b)
    if cancel_alpha_phs
        b_a_phase = fft(b) .* cis.(-angle.(fft(reverse(a)))) # todo: use sign
        b = ifft(b_a_phase)
    end
    return ab2rf(a, b)
end


"""
    b2a(b)

Convert an SLR beta parameter to its minimum-phase alpha parameter.

# Arguments
- `b::AbstractVector`: SLR beta parameter.

# Returns
- `Vector{ComplexF64}`: Minimum-phase SLR alpha parameter.
"""
function b2a(b::AbstractVector)
    n = length(b)
    npad = n * 16
    bcp = zeros(ComplexF64, npad)
    bcp[1:n] .= b
    bf = fft(bcp)
    bfmax = maximum(abs, bf)
    if bfmax ≥ 1
        @. bf /= (1e-7 + bfmax)
    end
    afa = mag2mp(@. sqrt(1 - abs2(bf)))
    a = fft(afa) ./ npad
    return reverse(a[1:n])
end


"""
    mag2mp(x)

Convert a magnitude spectrum `x` to its minimum-phase equivalent.

# Arguments
- `x::AbstractVector`: Magnitude spectrum.

# Returns
- `Vector{ComplexF64}`: Minimum-phase filter.
"""
function mag2mp(x::AbstractVector)
    n = length(x)
    xlf = fft(@. log(abs(x)))
    xlfp = copy(xlf)
    # Double positive frequencies, zero negative
    @. xlfp[2:n÷2] *= 2
    @. xlfp[n÷2 + 2:end] = 0
    return exp.(ifft(xlfp))
end


"""
    ab2rf(a, b)

Inverse SLR transform: alpha/beta polynomials to RF pulse.

# Arguments
- `a::AbstractVector`: SLR alpha parameter.
- `b::AbstractVector`: SLR beta parameter.

# Returns
- `Vector{ComplexF64}`: Designed RF pulse.
"""
function ab2rf(a::AbstractVector, b::AbstractVector)
    n = length(a)
    rf = zeros(ComplexF64, n)
    a, b = ComplexF64.(a), ComplexF64.(b)
    for ii in n:-1:1
        ratio = b[ii] / a[ii]
        cj = sqrt(1 / (1 + abs2(ratio)))
        sj = conj(cj * ratio)
        rf[ii] = 2 * atan(abs(sj), cj) * cis(angle(sj)) # todo sign
        if ii > 1
            at = @. cj * a + sj * b
            bt = @. -conj(sj) * a + cj * b
            a, b = at[2:ii], bt[1:ii-1]
        end
    end
    return rf
end


"""
    root_flip(b, d1, flip, tb; verbose=false)

Exhaustive root-flip pattern search for minimum peak B1.

# Arguments
- `b::AbstractVector`: SLR beta parameter.
- `d1::Real`: Passband ripple level.
- `flip::Real`: Target flip angle.
- `tb::Real`: pulse time bandwidth product.
- `verbose::Bool`: Print feedback on iterations.

# Returns
- `Tuple{Vector{ComplexF64}, Vector{Float64}}`: `(rf_out, b_out)`.

# References
- Sharma, A., Lustig, M. and Grissom, W. (2016). 'Root-flipped multiband refocusing pulses'.
  Magn. Reson. Med. 75(1), 227-237.
"""
function root_flip(b::AbstractVector, d1::Real, flip::Real, tb::Real; verbose::Bool=false)
    n = length(b)
    w = range(0, π, 512)
    b ./= maximum(abs, freqresp(PolynomialRatio(b, [1.0]), w))
    b .*= sin(flip / 2 + atan(d1 * 2) / 2)
    # Revert to standard descending orientation to match SigPy
    r = roots(Polynomial(reverse(b)))
    r = leja(r)
    candidates = @. (abs(1 - abs(r)) > 0.004) & (abs(angle(r)) < tb / n * π)
    n_cand = sum(candidates)
    rf_out, b_out, min_peak = b2rf(b), copy(b), maximum(abs, b2rf(b))
    for ii in 0:(2^n_cand - 1)
        # MSB first bit ordering to match Python's string-based search
        bits = reverse(digits(ii, base=2, pad=n_cand))
        r_f = copy(r)
        cand_indices = findall(candidates)
        for j in 1:n_cand
            if bits[j] == 1
                r_f[cand_indices[j]] = conj(1 / r_f[cand_indices[j]])
            end
        end
        b_t = reverse(coeffs(fromroots(r_f)))
        b_t ./= maximum(abs, freqresp(PolynomialRatio(b_t, [1.0]), w))
        b_t .*= sin(flip / 2 + atan(d1 * 2) / 2)
        rf_t = b2rf(b_t)
        if maximum(abs, rf_t) < min_peak
            rf_out, b_out, min_peak = rf_t, b_t, maximum(abs, rf_t)
        end
    end
    return rf_out, b_out
end


"""
    dz_recursive_rf(; n_seg, tb, n, ...)

Recursive SLR pulse design for variable flip angle or spin-echo sequences.

# Arguments
- `n_seg::Int`: Number of segments designed by recursion.
- `tb::Real`: Time-bandwidth product.
- `n::Int`: Pulse length.
- `se_seq::Bool`: Spin-echo sequence.
- `tb_ref::Real`: Time-bandwidth product of refocusing pulse.
- `z_pad_fact::Real`: Zero padding factor.
- `win_fact::Real`: Applied window factor.
- `cancel_alpha_phs::Bool`: Absorb alpha phase for flatter phase.
- `t1::Real`: T1 relaxation time.
- `tr_seg::Real`: Length of TR segment.
- `use_mz::Bool`: Design pulses accounting for actual Mz profile.
- `d1::Real`: Passband ripple level.
- `d2::Real`: Stopband ripple level.
- `d1se::Real`: Passband ripple level for spin-echo.
- `d2se::Real`: Stopband ripple level for spin-echo.

# Returns
- `Matrix{ComplexF64}` or `Tuple{Matrix{ComplexF64}, Vector{ComplexF64}}`: RF pulse(s).
"""
function dz_recursive_rf(;
    n_seg::Int, tb::Real, n::Int,
    se_seq::Bool = false, tb_ref::Real = 8, z_pad_fact::Real = 4.0, win_fact::Real = 1.75,
    cancel_alpha_phs::Bool = true, t1::Real = Inf, tr_seg::Real = 60.0, use_mz::Bool = true,
    d1::Real = 0.01, d2::Real = 0.01, d1se::Real = 0.01, d2se::Real = 0.01,
)
    # get refocusing pulse and its rotation parameters
    rf_ref = nothing
    bref = ones(ComplexF64, Int(z_pad_fact * n))
    bref_mag, aref_mag, flip_ref = zeros(length(bref)), zeros(length(bref)), 0.0

    if se_seq
        bsf, d1s, d2s = calc_ripples(:se, d1se, d2se)
        b_r = bsf .* dzls(n, tb_ref, d1s, d2s)
        pad = Int(z_pad_fact * n / 2 - n / 2)
        b_r_full = [zeros(pad); b_r; zeros(pad)]
        rf_ref = b2rf(b_r_full)
        bref = sp_fft(b_r_full) ./ maximum(abs, sp_fft(b_r_full))
        bref_mag, aref_mag = abs.(bref), @. sqrt(max(0, 1 - abs2(bref)))
        flip_ref = 2 * asin(bref_mag[Int(z_pad_fact * n / 2) + 1]) * 180 / π
    end

    # get flip angles (in degrees)
    flip = zeros(n_seg)
    flip[end] = 90.0
    for jj in n_seg-1:-1:1
        val = sin(flip[jj+1] * π / 180) * (se_seq ? cos(flip_ref * π / 180) : 1.0)
        flip[jj] = atan(val) * 180 / π
    end

    # design first RF pulse
    n_tot = Int(z_pad_fact * n)
    b = zeros(ComplexF64, n_tot, n_seg)
    idx = Int(z_pad_fact * n / 2 - n / 2) + 1 : Int(z_pad_fact * n / 2 + n / 2)
    b[idx, 1] .= dzls(n, tb, d1, d2)
    B = sp_fft(b[:, 1]) .* cis.(-2π / n_tot / 2 .* (-(n_tot/2) : (n_tot/2-1)))
    b[:, 1] .= sp_ifft(B ./ maximum(abs, B)) .* sin(flip[1] * π / 360)

    rf = zeros(ComplexF64, n_tot, n_seg)
    a = b2a(b[:, 1])
    if cancel_alpha_phs
        b_a_phase = fft(b[:, 1]) .* cis.(-angle.(fft(reverse(a)))) # todo: sign
        b[:, 1] .= ifft(b_a_phase)
    end
    rf[:, 1] .= ab2rf(a, b[:, 1])

    A, B = sp_fft(a), sp_fft(b[:, 1])
    window = ones(n_tot)
    if win_fact < z_pad_fact
        wl, np = Int((win_fact-1)*n), Int(n*z_pad_fact - win_fact*n)
        wc = blackman(wl)
        window = [zeros(np÷2); wc[1:wl÷2]; ones(n); wc[wl÷2+1:end]; zeros(np÷2)]
        length(window) != n_tot && (window = [window; zeros(n_tot - length(window))])
        b[:, 1] .*= window
        a = b2a(b[:, 1]) # Recalculate a after windowing
        rf[:, 1] .= ab2rf(a, b[:, 1])
        B, A = sp_fft(b[:, 1]), sp_fft(a)
    end

    mxy0 = se_seq ? @.(2 * A * conj(B) * bref^2) : @.(2 * conj(A) * B)
    mz = ones(ComplexF64, n_tot)
    for jj in 2:n_seg
        if se_seq
            @. mz *= (1 - 2 * (abs2(A * bref_mag) + abs2(aref_mag * B)))
        else
            @. mz = mz * (1 - 2 * abs2(B)) * exp(-tr_seg / t1) + (1 - exp(-tr_seg / t1))
        end

        if use_mz
            cq = @. -abs2(mxy0)
            bq = @. 4 * mz^2 * (se_seq ? bref_mag^4 : 1.0)
            aq = -bq
            delta = @. complex(bq^2 - 4 * aq * cq)
            # Use real part before max to avoid complex comparison
            bmag = @. sqrt(max(0.0, real((-bq + sqrt(delta)) / (2 * aq))))
            A = sp_fft(b2a(sp_ifft(bmag)))
            B = @. mxy0 / (2 * conj(A) * mz + 1e-16)
        else
            @. B *= sin(π/360 * flip[jj]) / sin(π/360 * flip[jj-1])
            A = sp_fft(b2a(sp_ifft(B)))
        end
        b[:, jj] .= sp_ifft(B) .* window
        rf[:, jj] .= b2rf(b[:, jj])
        B, A = sp_fft(b[:, jj]), sp_fft(b2a(b[:, jj]))
    end

    if win_fact < z_pad_fact
        pl = Int(win_fact * n)
        rf = rf[Int(n*z_pad_fact - win_fact*n)÷2 + 1 : Int(n*z_pad_fact - win_fact*n)÷2 + pl, :]
    end
    return se_seq ? (rf, rf_ref) : rf
end
