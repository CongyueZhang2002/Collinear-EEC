using PolyLog

if !isdefined(@__MODULE__, :CA)
    include("constants.jl")
end

# Timelike splitting kernels in the convention as = alpha_s / (4*pi).
# LO and NLO are exact. The NNLO expressions are the compact x-space
# parametrizations in Appendix B of arXiv:1107.2263.

function _pqq2_regular(x, nf, L0, L1)
    one_minus_x = 1 - x
    invx = inv(x)

    non_singlet = (
        1658.7 - 707.67 * L1 + 1327.5 * L0 - 56.907 * L0 * L1
        - 189.37 * L0^2 - 519.37 * L1 * L0^2
        - (352 / 9) * L0^3 + (128 / 81) * L0^4
        - 4249.4 * x - 559.1 * x * L0 * L1
        - 1075.3 * x^2 + 593.9 * x^3
        + nf * (
            (5120 / 81) * L1 - 198.1 + 466.29 * x
            + 181.18 * x^2 - 31.84 * x^3 - 39.113 * x * L0
            - L0 * L1 * (
                50.758 - 85.72 * x - 28.551 * L0 + 23.102 * x * L0
            )
            - 168.89 * L0 - (176 / 81) * L0^2 + (64 / 27) * L0^3
        )
        + nf^2 * (64 / 81) * (
            x * L0 * (1.5 * L0 + 5) / one_minus_x
            + 1
            + one_minus_x * (6 + 5.5 * L0 + 0.75 * L0^2)
        )
    )

    pure_singlet_nf = (
        -5.926 * L1^3 - 9.751 * L1^2 - 8.65 * L1
        - 106.65 - 848.97 * x + 368.79 * x^2 - 61.284 * x^3
        + 96.171 * L0 * L1 + 656.49 * L0 + 425.14 * L0^2
        + 47.322 * L0^3 + 9.072 * L0^4
        + 479.87 * invx + 324.07 * invx * L0
        - (128 / 9) * invx * L0^2 - (256 / 9) * invx * L0^3
    )

    pure_singlet_nf2 = (
        1.778 * L1^2 + 16.611 * L1 + 87.795
        - 57.688 * x - 41.827 * x^2 + 25.628 * x^3 - 7.9934 * x^4
        - 2.1031 * L0 * L1 + 26.294 * x * L0 - 7.8645 * x * L0^3
        + 57.713 * L0 + 9.1682 * L0^2 - 1.9 * L0^3
        + 0.019122 * L0^4 - (128 / 81) * invx
    )

    pure_singlet = one_minus_x * (
        nf * pure_singlet_nf + nf^2 * pure_singlet_nf2
    )

    return non_singlet + pure_singlet
end

function _pqg2(x, nf, L0, L1)
    invx = inv(x)

    nf1 = (
        (100 / 27) * L1^4 + (350 / 9) * L1^3
        + 263.07 * L1^2 + 693.84 * L1 + 603.71
        - 882.48 * x + 4723.2 * x^2 - 4745.8 * x^3 - 175.28 * x^4
        - L0 * L1 * (1809.4 + 107.59 * x)
        - 885.5 * x * L0^4
        + 1864 * L0 + 1512 * L0^2 + 361.28 * L0^3 + 42.328 * L0^4
        + 1141.7 * invx + 675.83 * invx * L0
        - 64 * invx * (L0^2 + L0^3)
    )

    nf2 = (
        -(100 / 27) * L1^3 - 35.446 * L1^2 - 103.609 * L1
        - 113.81 + 341.26 * x - 853.35 * x^2 + 492.1 * x^3
        + 14.803 * x^4
        + L0 * L1 * (966.96 - 1.593 * L1 - 709.1 * x)
        - 333.8 * x * L0^3
        + 619.75 * L0 + 255.62 * L0^2 + 21.569 * L0^3
        - 2.8986 * invx - 3.1752 * invx * L0
        - (32 / 27) * invx * L0^2
    )

    nf3 = (4 / 9) * (
        4 + 6 * (L0 + L1)
        + (1 - 2 * x + 2 * x^2) * (
            3.8696 + 4 * (L0 + L1) + 3 * (L0 + L1)^2
        )
    )

    return nf * nf1 + nf^2 * nf2 + nf^3 * nf3
end

function _pgq2(x, nf, L0, L1)
    invx = inv(x)

    nf0 = (
        (400 / 81) * L1^4 + (520 / 27) * L1^3
        - 220.13 * L1^2 - 152.6 * L1 + 272.85
        - 7188.7 * x + 5693.2 * x^2 + 146.98 * x^3 + 128.19 * x^4
        - L0 * L1 * (1300.6 + 71.23 * L1)
        + 543.8 * x * L0^3
        + 4.4136 * L0 - 0.71252 * L0^2 - 126.38 * L0^3
        - 30.061 * L0^4
        + 5803.7 * invx + 4776.5 * invx * L0
        + 1001.89 * invx * L0^2 + (3712 / 3) * invx * L0^3
        + 256 * invx * L0^4
    )

    nf1 = (
        (80 / 81) * L1^3 + (1040 / 81) * L1^2 - 16.914 * L1
        - 871.3 + 790.13 * x - 241.23 * x^2 + 43.252 * x^3
        - 4.3465 * x * L0^3 + 55.048 * L0 * L1
        - 492 * L0 - 343.1 * L0^2 - 48.6 * L0^3
        + 6.0041 * invx + 141.93 * invx * L0
        + (2912 / 27) * invx * L0^2 + (1280 / 81) * invx * L0^3
    )

    return nf0 + nf * nf1
end

function _pgg2_regular(x, nf, L0, L1)
    invx = inv(x)

    nf0 = (
        -3590.1 * L1 - 28489 + 7469 * x + 30421 * x^2
        - 53017 * x^3 + 19556 * x^4
        - L0 * L1 * (186.4 + 21328 * L0)
        + 12258 * L0 + 13528 * L0^2 + 3281.7 * L0^3
        + 191.99 * L0^4 + 5685.8 * x * L0^3
        + 14214.4 * invx + 10233 * invx * L0
        + 3651.1 * invx * L0^2 + 3168 * invx * L0^3
        + 576 * invx * L0^4
    )

    nf1 = (
        319.97 * L1 + 248.95 + 260.6 * x + 272.79 * x^2
        + 2133.2 * x^3 - 926.87 * x^4
        + L0 * L1 * (1266.5 - 29.709 * L0 + 87.771 * L1)
        + 4.9934 * L0 + 482.94 * L0^2 + 155.1 * L0^3
        + 18.085 * L0^4 + 485.18 * x * L0^3
        - 804.13 * invx - 5.47 * invx * L0
        + (2368 / 9) * invx * L0^2 + (448 / 9) * invx * L0^3
    )

    nf2 = (
        -77.19 + 153.27 * x - 106.03 * x^2 + 11.995 * x^3
        - L0 * L1 * (115.01 - 96.522 * x + 62.908 * L0)
        - 69.712 * L0 - 44.8 * L0^2 - 5.037 * L0^3
        + (472 / 243) * invx + (368 / 81) * invx * L0
        + (32 / 27) * invx * L0^2
    )

    return nf0 + nf * nf1 + nf^2 * nf2
end

"""
    timelike_splittings_func(; x, type, nf)

Return one timelike splitting-kernel coefficient at the requested loop order.
Diagonal channels return `(regular, delta, D0)` at LO/NNLO and
`(regular, delta, D0, D1)` at NLO. Off-diagonal channels return a scalar.
"""
function timelike_splittings_func(; x::Float64, type::String, nf::Int64)
    if !(0.0 < x < 1.0)
        throw(DomainError(x, "x must lie strictly between 0 and 1."))
    end
    if nf < 0
        throw(DomainError(nf, "nf must be nonnegative."))
    end

    x2 = x^2
    x3 = x^3
    x4 = x^4
    invx = inv(x)
    inv1px = inv(1 + x)

    L0 = log(x)
    L1 = log1p(-x)

    if type == "Pqq0"
        return (-(8 / 3) * (1 + x), 4.0, 16 / 3)
    elseif type == "Pqg0"
        return nf * (2 - 4 * x + 4 * x2)
    elseif type == "Pgq0"
        return (4 / 3) * (-4 + 4 * invx + 2 * x)
    elseif type == "Pgg0"
        return (
            12 * (-2 + invx + x - x2),
            11 - (2 / 3) * nf,
            12.0,
        )
    end

    if type == "Pqq2"
        regular = _pqq2_regular(x, nf, L0, L1)
        delta = (
            1295.625 - 173.935 * nf
            - (64 / 81) * nf^2 * (51 / 16 + 3 * z3 - 5 * z2)
        )
        D0 = 1174.898 - 183.187 * nf - (64 / 81) * nf^2
        return (regular, delta, D0)
    elseif type == "Pqg2"
        return _pqg2(x, nf, L0, L1)
    elseif type == "Pgq2"
        return _pgq2(x, nf, L0, L1)
    elseif type == "Pgg2"
        regular = _pgg2_regular(x, nf, L0, L1)
        delta = 4425.451 - 528.719 * nf + 6.4628 * nf^2
        D0 = 2643.521 - 412.172 * nf - (16 / 9) * nf^2
        return (regular, delta, D0)
    end

    H0 = L0
    H1 = -L1
    H2 = reli(2, x)
    Hm10 = L0 * log1p(x) + reli(2, -x)
    H00 = 0.5 * L0^2
    H10 = reli(2, 1 - x) - z2
    H11 = 0.5 * L1^2

    if type == "Pqq1"
        regular = (8 / 27) * (
            27 - 146 * nf + 8 * pi^2 - 40 * nf * invx
            - 429 * x + 94 * nf * x + 10 * pi^2 * x + 112 * nf * x2
            + 2 * pi^2 * inv1px
            - 3 * (57 + 41 * x + 4 * nf * (7 + 13 * x + 4 * x2)) * H0
            + 12 * (6 + 5 * x - inv1px + 3 * nf * (1 + x)) * H00
            + 48 * (1 + x) * (H10 + H2)
            + (-12 + 12 * x + 24 * inv1px) * Hm10
        )
        delta = (-2 / 27) * (
            -189 - 84 * pi^2 + nf * (6 + 8 * pi^2) + 72 * z3
        )
        D0 = (-16 / 27) * (
            -201 + 10 * nf + 9 * pi^2
            + 3 * (-45 + 2 * nf) * L0 + 21 * L0^2
        )
        D1 = (256 / 9) * L0
        return (regular, delta, D0, D1)
    elseif type == "Pqg1"
        return (-4 * nf / (9 * x)) * (
            60
            + x * (
                33 + 3 * (49 - 138 * x) * x
                + 2 * nf * (5 + 4 * (x - 1) * x)
                + 2 * pi^2 * (2 + x * (5 + 4 * x))
            )
            + 6 * x * (
                11 + nf + 2 * nf * (x - 1) * x + x * (47 + 2 * x)
            ) * H0
            + 6 * x * (-7 - 58 * x + 8 * x2) * H00
            + 3 * x * (
                (21 + 22 * (x - 1) * x - 2 * nf * (1 + 2 * (x - 1) * x)) * H1
                - 2 * (1 + 2 * (x - 1) * x) * (6 * H10 + 5 * H11 - 14 * H2)
                + 18 * (1 + 2 * x * (1 + x)) * Hm10
            )
        )
    elseif type == "Pgq1"
        return (16 / (9 * x)) * (
            17 + x * (43 + 6 * pi^2 + (9 - 44 * x) * x)
            + (-54 + x * (40 + x * (83 + 24 * x))) * H0
            - 2 * (36 + x * (14 + 29 * x)) * H00
            - 76 * H10 - 20 * H11 + 4 * H2 + 36 * Hm10
            + 2 * x * (
                5 * x * H1
                - (x - 2) * (19 * H10 + 5 * H11 - H2)
                + 9 * (2 + x) * Hm10
            )
        )
    elseif type == "Pgg1"
        regular = (
            -4 * nf * (1 + x) * (23 + x * (-189 + x * (-45 + 121 * x)))
            + 54 * x * (
                -(1 + x) * (25 + 109 * x)
                + 6 * pi^2 * (3 + 2 * x * (2 + x + x2))
            )
            + 12 * (1 + x) * (
                -27 * (22 + x * (33 + x * (3 + 22 * x)))
                + 2 * nf * (-2 + x * (57 + x * (15 + 34 * x)))
            ) * H0
            - 3888 * H10 - 3888 * H2 + 3888 * Hm10
            + 72 * (
                (
                    -108 + x * (
                        -81 + 54 * (x - 4) * x * (1 + x)
                        + 4 * nf * (1 + x)^2
                    )
                ) * H00
                + 54 * x * (
                    (1 + x + x3) * (H10 + H2)
                    + (1 + x) * (2 + x + x2) * Hm10
                )
            )
        ) / (27 * x * (1 + x))
        delta = 96 - (32 / 3) * nf + 108 * z3
        D0 = (
            8 * (33 - 2 * nf) * L0
            - (4 / 3) * (-201 + 10 * nf + 9 * pi^2 + 81 * L0^2)
        )
        D1 = 144 * L0
        return (regular, delta, D0, D1)
    end

    throw(ArgumentError("unknown splitting-kernel type: $type"))
end

"""
    timelike_splitting_convolution_func(; y, as, order, nf)

Return the timelike singlet kernel matrix components through `order`, including
the appropriate powers of `as = alpha_s / (4*pi)`.
"""
function timelike_splitting_convolution_func(;
    y::Float64,
    as::Float64,
    order::Int64,
    nf::Int64,
)
    if !(0.0 < y < 1.0)
        throw(DomainError(y, "y must lie strictly between 0 and 1."))
    end
    if !isfinite(as) || as < 0.0
        throw(DomainError(as, "as must be finite and nonnegative."))
    end
    if !(0 <= order <= 2)
        throw(ArgumentError("order must be 0, 1, or 2."))
    end
    if nf < 0
        throw(DomainError(nf, "nf must be nonnegative."))
    end

    Pqq0Reg, Pqq0Delta, Pqq0D0 = timelike_splittings_func(
        x = y, type = "Pqq0", nf = nf,
    )
    Pqg0 = timelike_splittings_func(x = y, type = "Pqg0", nf = nf)
    Pgq0 = timelike_splittings_func(x = y, type = "Pgq0", nf = nf)
    Pgg0Reg, Pgg0Delta, Pgg0D0 = timelike_splittings_func(
        x = y, type = "Pgg0", nf = nf,
    )

    PqqReg = as * Pqq0Reg
    PqqDelta_at_1 = as * Pqq0Delta
    PqqD0 = as * Pqq0D0
    PqqD0_at_1 = as * Pqq0D0
    PqqD1 = 0.0
    PqqD1_at_1 = 0.0
    Pqg = as * Pqg0
    Pgq = as * Pgq0

    PggReg = as * Pgg0Reg
    PggDelta_at_1 = as * Pgg0Delta
    PggD0 = as * Pgg0D0
    PggD0_at_1 = as * Pgg0D0
    PggD1 = 0.0
    PggD1_at_1 = 0.0

    if order >= 1
        as1 = as^2
        Pqq1Reg, Pqq1Delta, Pqq1D0, Pqq1D1 = timelike_splittings_func(
            x = y, type = "Pqq1", nf = nf,
        )
        Pqg1 = timelike_splittings_func(x = y, type = "Pqg1", nf = nf)
        Pgq1 = timelike_splittings_func(x = y, type = "Pgq1", nf = nf)
        Pgg1Reg, Pgg1Delta, Pgg1D0, Pgg1D1 = timelike_splittings_func(
            x = y, type = "Pgg1", nf = nf,
        )

        PqqReg += as1 * Pqq1Reg
        PqqDelta_at_1 += as1 * Pqq1Delta
        PqqD0 += as1 * Pqq1D0
        PqqD0_at_1 += as1 * (
            (-16 / 27) * (-201 + 10 * nf + 9 * pi^2)
        )
        PqqD1 += as1 * Pqq1D1
        Pqg += as1 * Pqg1
        Pgq += as1 * Pgq1

        PggReg += as1 * Pgg1Reg
        PggDelta_at_1 += as1 * Pgg1Delta
        PggD0 += as1 * Pgg1D0
        PggD0_at_1 += as1 * (
            (-4 / 3) * (-201 + 10 * nf + 9 * pi^2)
        )
        PggD1 += as1 * Pgg1D1
    end

    if order == 2
        as2 = as^3
        Pqq2Reg, Pqq2Delta, Pqq2D0 = timelike_splittings_func(
            x = y, type = "Pqq2", nf = nf,
        )
        Pqg2 = timelike_splittings_func(x = y, type = "Pqg2", nf = nf)
        Pgq2 = timelike_splittings_func(x = y, type = "Pgq2", nf = nf)
        Pgg2Reg, Pgg2Delta, Pgg2D0 = timelike_splittings_func(
            x = y, type = "Pgg2", nf = nf,
        )

        PqqReg += as2 * Pqq2Reg
        PqqDelta_at_1 += as2 * Pqq2Delta
        PqqD0 += as2 * Pqq2D0
        PqqD0_at_1 += as2 * Pqq2D0
        Pqg += as2 * Pqg2
        Pgq += as2 * Pgq2

        PggReg += as2 * Pgg2Reg
        PggDelta_at_1 += as2 * Pgg2Delta
        PggD0 += as2 * Pgg2D0
        PggD0_at_1 += as2 * Pgg2D0
    end

    return (
        PqqReg = PqqReg,
        PqqDelta_at_1 = PqqDelta_at_1,
        PqqD0 = PqqD0,
        PqqD0_at_1 = PqqD0_at_1,
        PqqD1 = PqqD1,
        PqqD1_at_1 = PqqD1_at_1,
        Pqg = Pqg,
        Pgq = Pgq,
        PggReg = PggReg,
        PggDelta_at_1 = PggDelta_at_1,
        PggD0 = PggD0,
        PggD0_at_1 = PggD0_at_1,
        PggD1 = PggD1,
        PggD1_at_1 = PggD1_at_1,
    )
end

"""
    timelike_splitting_coefficient_func(; y, loop_order, nf)

Return `P^(loop_order)(y)` without a power of `as`.
"""
function timelike_splitting_coefficient_func(;
    y::Float64,
    loop_order::Int64,
    nf::Int64,
)
    if !(0 <= loop_order <= 2)
        throw(ArgumentError("loop_order must be 0, 1, or 2."))
    end

    current = timelike_splitting_convolution_func(
        y = y,
        as = 1.0,
        order = loop_order,
        nf = nf,
    )
    loop_order == 0 && return current

    previous = timelike_splitting_convolution_func(
        y = y,
        as = 1.0,
        order = loop_order - 1,
        nf = nf,
    )
    names = propertynames(current)
    return NamedTuple{names}(
        ntuple(
            index -> begin
                name = names[index]
                getproperty(current, name) - getproperty(previous, name)
            end,
            length(names),
        ),
    )
end
