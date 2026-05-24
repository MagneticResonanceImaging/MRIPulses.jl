# util.jl


"""
    dinf(d1 = 0.01, d2 = 0.01)

Calculate `D_∞` for a linear phase filter.

# Arguments
- `d1::Real`: passband ripple level in M0⁻¹.
- `d2::Real`: stopband ripple level in M0⁻¹.

# Returns
- `Real`: D infinity.

# References
Pauly, Le Roux, Nishimura, Macovski:
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm.
IEEE Tr Medical Imaging 1991; 10(1):53-65.
https://doi.org/10.1109/42.75611
"""
function dinf(d1::Real = 0.01, d2::Real = 0.01)
    a1 = 5.309e-3
    a2 = 7.114e-2
    a3 = -4.761e-1
    a4 = -2.66e-3
    a5 = -5.941e-1
    a6 = -4.278e-1

    l10d1 = log10(d1)
    l10d2 = log10(d2)

    d = (a1 * l10d1^2 + a2 * l10d1 + a3) * l10d2 +
        (a4 * l10d1^2 + a5 * l10d1 + a6)

    return d
end
