# Beta-function and cusp coefficients adapted from the WSL TMD fitter
# Core/anomalous dims.jl. The beta coefficients are evaluated as functions of
# nf so the numerical coupling is not tied to a particular flavor count.

if !isdefined(@__MODULE__, :CA)
    include("constants.jl")
end

function _check_nf(nf::Integer)
    nf >= 0 || throw(DomainError(nf, "nf must be nonnegative."))
    return Int(nf)
end

β0_func(nf::Integer) = 11.0 - (2.0 / 3.0) * _check_nf(nf)

β1_func(nf::Integer) = 102.0 - (38.0 / 3.0) * _check_nf(nf)

function β2_func(nf::Integer)
    n = _check_nf(nf)
    return 2857.0 / 2.0 - (5033.0 / 18.0) * n + (325.0 / 54.0) * n^2
end

function β3_func(nf::Integer)
    n = _check_nf(nf)
    return (
        149753.0 / 6.0 + 3564.0 * z3
        - (1078361.0 / 162.0 + (6508.0 / 27.0) * z3) * n
        + (50065.0 / 162.0 + (6472.0 / 81.0) * z3) * n^2
        + (1093.0 / 729.0) * n^3
    )
end

function β4_func(nf::Integer)
    n = _check_nf(nf)
    return (
        8157455.0 / 16.0 + (621885.0 / 2.0) * z3
        - (88209.0 / 2.0) * z4 - 288090.0 * z5
        + n * (
            -336460813.0 / 1944.0 - (4811164.0 / 81.0) * z3
            + (33935.0 / 6.0) * z4 + (1358995.0 / 27.0) * z5
        )
        + n^2 * (
            25960913.0 / 1944.0 + (698531.0 / 81.0) * z3
            - (10526.0 / 9.0) * z4 - (381760.0 / 81.0) * z5
        )
        + n^3 * (
            -630559.0 / 5832.0 - (48722.0 / 243.0) * z3
            + (1618.0 / 27.0) * z4 + (460.0 / 9.0) * z5
        )
        + n^4 * (1205.0 / 2916.0 - (152.0 / 81.0) * z3)
    )
end

Γ0_func(nf::Integer) = (_check_nf(nf); 4.0 * CF)

function Γ1_func(nf::Integer)
    n = _check_nf(nf)
    return 4.0 * CF * (
        (67.0 / 9.0 - pi^2 / 3.0) * CA
        - (20.0 / 9.0) * TF * n
    )
end

function Γ2_func(nf::Integer)
    n = _check_nf(nf)
    return 4.0 * CF * (
        CA^2 * (
            245.0 / 6.0 - 134.0 * pi^2 / 27.0
            + 11.0 * pi^4 / 45.0 + (22.0 / 3.0) * z3
        )
        + CA * TF * n * (
            -418.0 / 27.0 + 40.0 * pi^2 / 27.0
            - (56.0 / 3.0) * z3
        )
        + CF * TF * n * (-55.0 / 3.0 + 16.0 * z3)
        - (16.0 / 27.0) * TF^2 * n^2
    )
end

# The copied fitter stores these higher-order values numerically.
function Γ3_func(nf::Integer)
    nf == 3 && return 7.035152974055e3
    nf == 4 && return 3.353354170411e3
    nf == 5 && return 1.412460850831e2
    throw(ArgumentError("Γ3_func currently has tabulated values only for nf = 3, 4, or 5."))
end

function Γ4_func(nf::Integer)
    nf in 3:5 ||
        throw(ArgumentError("Γ4_func currently has an estimate only for nf = 3, 4, or 5."))
    return 0.21
end

"""
    Γ_func(; mu, order, nf_scheme=:VFNS, ...)

Evaluate the copied TMD fitter cusp series using the standalone numerical
coupling. The fitter convention uses `order = 1` through `4` for truncation
after the two- through five-loop terms, respectively.
"""
function Γ_func(;
    mu::Real,
    order::Integer,
    nf_scheme::Symbol = :VFNS,
    alpha_s_ref::Real = DEFAULT_ALPHA_S_MZ,
    alpha_s_mu_ref::Real = MZ,
    alpha_s_order::Integer = 4,
)
    1 <= order <= 4 || throw(ArgumentError("order must be between 1 and 4."))
    nf = nf_func(mu; scheme = nf_scheme)
    alpha_s = alpha_s_func_numerical(
        mu_f = mu,
        mu_ref = alpha_s_mu_ref,
        alpha_s_ref = alpha_s_ref,
        order = alpha_s_order,
        nf_scheme = nf_scheme,
    )
    a_s = alpha_s / (4.0 * pi)
    value = a_s * Γ0_func(nf) + a_s^2 * Γ1_func(nf)
    order >= 2 && (value += a_s^3 * Γ2_func(nf))
    order >= 3 && (value += a_s^4 * Γ3_func(nf))
    order >= 4 && (value += a_s^5 * Γ4_func(nf))
    return value
end
