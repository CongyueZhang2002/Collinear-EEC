# Reference: Appendix A, Eqs. (A4)-(A12), arXiv:1905.01310.
# The paper defines gamma_T(N=3) = -integral_0^1 dx x^2 P_T(x).
using QuadGK

const PAPER_ZETA2 = pi^2 / 6.0
const PAPER_ZETA3 = 1.2020569031595943
const PAPER_ZETA4 = pi^4 / 90.0
const PAPER_ZETA5 = 1.03692775514337

if !isdefined(@__MODULE__, :timelike_splitting_coefficient_func)
    include(joinpath(@__DIR__, "..", "core", "timelike_splitting.jl"))
end

function loop_splitting_components(
    y::Float64,
    loop_order::Int;
    nf::Int,
)
    return timelike_splitting_coefficient_func(
        y = y,
        loop_order = Int64(loop_order),
        nf = Int64(nf),
    )
end

function integrate_unit_interval(integrand)
    return quadgk(
        integrand,
        0.0,
        0.5,
        0.9,
        0.99,
        0.9999,
        1.0,
        rtol = 2.0e-10,
    )[1]
end

function splitting_moment(
    loop_order::Int;
    mellin_N::Int,
    log_power::Int = 0,
    nf::Int,
)

    mellin_N >= 1 || throw(ArgumentError("mellin_N must be positive."))
    log_power >= 0 || throw(ArgumentError("log_power must be nonnegative."))
    endpoint = loop_splitting_components(0.5, loop_order; nf = nf)
    weight_at_1 = log_power == 0 ? 1.0 : 0.0

    function weight(y)
        logarithm = log_power == 0 ? 1.0 : log(y)^log_power
        return y^(mellin_N - 1) * logarithm
    end

    function diagonal_moment(
        regular_name::Symbol,
        delta_name::Symbol,
        d0_name::Symbol,
        d0_at_1_name::Symbol,
        d1_name::Symbol,
        d1_at_1_name::Symbol,
    )
        function integrand(y)
            components = loop_splitting_components(y, loop_order; nf = nf)
            weighted_value = weight(y)
            return (
                weighted_value * getproperty(components, regular_name) +
                (
                    weighted_value * getproperty(components, d0_name) -
                    weight_at_1 * getproperty(endpoint, d0_at_1_name)
                ) / (1.0 - y) +
                (
                    weighted_value * getproperty(components, d1_name) -
                    weight_at_1 * getproperty(endpoint, d1_at_1_name)
                ) * log1p(-y) / (1.0 - y)
            )
        end

        delta = weight_at_1 * getproperty(endpoint, delta_name)
        return -(delta + integrate_unit_interval(integrand))
    end

    qq = diagonal_moment(
        :PqqReg,
        :PqqDelta_at_1,
        :PqqD0,
        :PqqD0_at_1,
        :PqqD1,
        :PqqD1_at_1,
    )
    gg = diagonal_moment(
        :PggReg,
        :PggDelta_at_1,
        :PggD0,
        :PggD0_at_1,
        :PggD1,
        :PggD1_at_1,
    )
    qg = -integrate_unit_interval(
        y -> weight(y) * loop_splitting_components(y, loop_order; nf = nf).Pqg,
    )
    gq = -integrate_unit_interval(
        y -> weight(y) * loop_splitting_components(y, loop_order; nf = nf).Pgq,
    )

    return (qq = qq, qg = qg, gq = gq, gg = gg)
end

function splitting_moment_N3(
    loop_order::Int;
    log_power::Int = 0,
    nf::Int,
)
    return splitting_moment(
        loop_order,
        mellin_N = 3,
        log_power = log_power,
        nf = nf,
    )
end

function paper_N3_references(nf_value::Real)
    ca = 3.0
    cf = 4.0 / 3.0
    nf = Float64(nf_value)
    z2p = PAPER_ZETA2
    z3p = PAPER_ZETA3
    z4p = PAPER_ZETA4
    z5p = PAPER_ZETA5

    lo = (
        qq = (25.0 / 6.0) * cf,
        qg = -(7.0 / 15.0) * nf,
        gq = -(7.0 / 6.0) * cf,
        gg = (14.0 / 5.0) * ca + (2.0 / 3.0) * nf,
    )
    nlo = (
        qq = (
            (-16.0 * z3p + 24.0 * z2p - 1693.0 / 48.0) * cf^2 +
            (8.0 * z3p - (86.0 / 3.0) * z2p + 459.0 / 8.0) * ca * cf -
            (5453.0 / 1800.0) * cf * nf
        ),
        qg = (
            ((28.0 / 15.0) * z2p + 619.0 / 2700.0) * ca * nf -
            (833.0 / 216.0) * cf * nf -
            (4.0 / 25.0) * nf^2
        ),
        gq = (
            ((28.0 / 3.0) * z2p - 2977.0 / 432.0) * cf^2 +
            (-(14.0 / 3.0) * z2p - 39451.0 / 5400.0) * ca * cf
        ),
        gg = (
            (-8.0 * z3p + (52.0 / 15.0) * z2p + 2158.0 / 675.0) * ca^2 +
            (-(16.0 / 3.0) * z2p + 3803.0 / 1350.0) * ca * nf +
            (12839.0 / 5400.0) * cf * nf
        ),
    )
    nnlo = (
        qq = (
            (
                112.0 * z5p + 48.0 * z2p * z3p - (2083.0 / 3.0) * z4p +
                (16153.0 / 18.0) * z3p - (13105.0 / 72.0) * z2p -
                3049531.0 / 31104.0
            ) * cf * ca^2 +
            (
                -432.0 * z5p - 208.0 * z2p * z3p + (8252.0 / 3.0) * z4p -
                (19424.0 / 9.0) * z3p - (16709.0 / 27.0) * z2p +
                20329835.0 / 15552.0
            ) * cf^2 * ca +
            (
                416.0 * z5p + 224.0 * z2p * z3p - (6172.0 / 3.0) * z4p +
                (10942.0 / 9.0) * z3p + (11797.0 / 18.0) * z2p -
                17471825.0 / 15552.0
            ) * cf^3 +
            (
                (68.0 / 3.0) * z4p - (5803.0 / 45.0) * z3p +
                (146971.0 / 2700.0) * z2p - 25234031.0 / 1944000.0
            ) * ca * cf * nf +
            (
                -(136.0 / 3.0) * z4p + (8176.0 / 45.0) * z3p -
                (9767.0 / 225.0) * z2p - 4100189.0 / 64800.0
            ) * cf^2 * nf -
            (105799.0 / 162000.0) * cf * nf^2
        ),
        qg = (
            (
                -(252.0 / 5.0) * z4p + (343.0 / 45.0) * z3p +
                (239959.0 / 13500.0) * z2p - 1795237.0 / 1944000.0
            ) * ca^2 * nf +
            (
                -(42.0 / 5.0) * z4p + (6208.0 / 75.0) * z3p +
                (34127.0 / 1350.0) * z2p - 3607891.0 / 38880.0
            ) * ca * cf * nf +
            (
                (448.0 / 15.0) * z4p - (26102.0 / 225.0) * z3p -
                (2042.0 / 225.0) * z2p + 9397651.0 / 97200.0
            ) * cf^2 * nf +
            (
                -(28.0 / 9.0) * z3p - (554.0 / 135.0) * z2p +
                1215691.0 / 121500.0
            ) * ca * nf^2 +
            ((2738.0 / 675.0) * z2p - 10657.0 / 4050.0) * cf * nf^2 -
            (172.0 / 1125.0) * nf^3
        ),
        gq = (
            (
                (196.0 / 3.0) * z4p - (2791.0 / 90.0) * z3p -
                (50593.0 / 600.0) * z2p - 17093053.0 / 777600.0
            ) * cf * ca^2 +
            (
                (511.0 / 3.0) * z4p - (3029.0 / 9.0) * z3p +
                (123773.0 / 900.0) * z2p + 63294389.0 / 388800.0
            ) * cf^2 * ca +
            (
                -308.0 * z4p + (2533.0 / 9.0) * z3p +
                (3193.0 / 54.0) * z2p - 647639.0 / 3888.0
            ) * cf^3 +
            (
                (182.0 / 9.0) * z3p - (73.0 / 27.0) * z2p +
                246767.0 / 60750.0
            ) * ca * cf * nf +
            (
                -(28.0 / 9.0) * z3p + (4.0 / 9.0) * z2p -
                419593.0 / 81000.0
            ) * cf^2 * nf
        ),
        gg = (
            (
                96.0 * z5p + 64.0 * z2p * z3p - (2566.0 / 15.0) * z4p -
                (23702.0 / 225.0) * z3p + (66358.0 / 1125.0) * z2p -
                5819653.0 / 486000.0
            ) * ca^3 +
            (
                104.0 * z4p + (239.0 / 9.0) * z3p -
                (51269.0 / 540.0) * z2p - 12230737.0 / 1944000.0
            ) * ca^2 * nf +
            (
                (282.0 / 5.0) * z3p - (16291.0 / 675.0) * z2p -
                1700563.0 / 108000.0
            ) * ca * cf * nf +
            (
                -(28.0 / 9.0) * z3p + (2411.0 / 675.0) * z2p +
                219077.0 / 194400.0
            ) * cf^2 * nf +
            (
                -(64.0 / 9.0) * z3p + (160.0 / 27.0) * z2p -
                18269.0 / 10125.0
            ) * ca * nf^2 +
            (
                -(196.0 / 135.0) * z2p - 2611.0 / 162000.0
            ) * cf * nf^2
        ),
    )

    return (lo, nlo, nnlo)
end

function paper_log_moment_references(nf_value::Real)
    ca = 3.0
    cf = 4.0 / 3.0
    nf = Float64(nf_value)
    z2p = PAPER_ZETA2
    z3p = PAPER_ZETA3
    z4p = PAPER_ZETA4

    lo_first = (
        qq = (4.0 * z2p - 385.0 / 72.0) * cf,
        qg = (119.0 / 900.0) * nf,
        gq = (49.0 / 72.0) * cf,
        gg = (4.0 * z2p - 4319.0 / 900.0) * ca,
    )
    lo_second = (
        qq = (-8.0 * z3p + 3979.0 / 432.0) * cf,
        qg = -(2353.0 / 27000.0) * nf,
        gq = -(331.0 / 432.0) * cf,
        gg = (-8.0 * z3p + 230353.0 / 27000.0) * ca,
    )
    nlo_first = (
        qq = (
            (
                -56.0 * z4p - (158.0 / 3.0) * z3p +
                (385.0 / 18.0) * z2p + 152863.0 / 1728.0
            ) * cf^2 +
            (
                -12.0 * z4p + (41.0 / 3.0) * z3p +
                (307.0 / 6.0) * z2p - 35785.0 / 432.0
            ) * cf * ca +
            (
                (16.0 / 3.0) * z3p - (40.0 / 9.0) * z2p -
                101923.0 / 108000.0
            ) * cf * nf
        ),
        qg = (
            (
                (42.0 / 5.0) * z3p - (92.0 / 75.0) * z2p -
                1460321.0 / 162000.0
            ) * ca * nf +
            (
                -(28.0 / 3.0) * z3p + (178.0 / 225.0) * z2p +
                46663.0 / 4320.0
            ) * cf * nf +
            (
                -(28.0 / 45.0) * z2p + 18451.0 / 20250.0
            ) * nf^2
        ),
        gq = (
            (
                -(49.0 / 3.0) * z3p + (59.0 / 6.0) * z2p +
                956963.0 / 108000.0
            ) * cf * ca +
            (
                14.0 * z3p - (275.0 / 18.0) * z2p +
                8053.0 / 1728.0
            ) * cf^2
        ),
        gg = (
            (
                -68.0 * z4p - (686.0 / 15.0) * z3p +
                (15338.0 / 225.0) * z2p + 3642257.0 / 162000.0
            ) * ca^2 +
            (
                (32.0 / 3.0) * z3p - (40.0 / 9.0) * z2p -
                137323.0 / 20250.0
            ) * ca * nf -
            (58247.0 / 108000.0) * cf * nf
        ),
    )

    return (lo_first = lo_first, lo_second = lo_second, nlo_first = nlo_first)
end

function test_moment_channels(actual, expected; rtol, atol)
    for channel in propertynames(expected)
        @test getproperty(actual, channel) ≈ getproperty(expected, channel) rtol = rtol atol = atol
    end
end

@testset "Published N=3 splitting moments" begin
    nf = 5
    references = paper_N3_references(nf)
    test_moment_channels(
        splitting_moment_N3(0; nf = nf),
        references[1],
        rtol = 2.0e-9,
        atol = 2.0e-9,
    )
    test_moment_channels(
        splitting_moment_N3(1; nf = nf),
        references[2],
        rtol = 2.0e-8,
        atol = 2.0e-7,
    )
    test_moment_channels(
        splitting_moment_N3(2; nf = nf),
        references[3],
        rtol = 1.0e-3,
        atol = 1.0e-2,
    )
end

@testset "Published logarithmic splitting moments" begin
    nf = 5
    references = paper_log_moment_references(nf)
    test_moment_channels(
        splitting_moment_N3(0; log_power = 1, nf = nf),
        references.lo_first,
        rtol = 2.0e-9,
        atol = 2.0e-9,
    )
    test_moment_channels(
        splitting_moment_N3(0; log_power = 2, nf = nf),
        references.lo_second,
        rtol = 2.0e-9,
        atol = 2.0e-9,
    )
    test_moment_channels(
        splitting_moment_N3(1; log_power = 1, nf = nf),
        references.nlo_first,
        rtol = 2.0e-8,
        atol = 2.0e-7,
    )
end

@testset "nf-dependent timelike splitting moments" begin
    for nf in 3:6
        n3_references = paper_N3_references(nf)
        test_moment_channels(
            splitting_moment_N3(0; nf = nf),
            n3_references[1],
            rtol = 2.0e-9,
            atol = 2.0e-9,
        )
        test_moment_channels(
            splitting_moment_N3(1; nf = nf),
            n3_references[2],
            rtol = 2.0e-8,
            atol = 2.0e-7,
        )
        # The source-faithful NNLO kernels are compact numerical parametrizations.
        test_moment_channels(
            splitting_moment_N3(2; nf = nf),
            n3_references[3],
            rtol = 2.0e-4,
            atol = 5.0e-2,
        )

        logarithmic_references = paper_log_moment_references(nf)
        test_moment_channels(
            splitting_moment_N3(0; log_power = 1, nf = nf),
            logarithmic_references.lo_first,
            rtol = 2.0e-9,
            atol = 2.0e-9,
        )
        test_moment_channels(
            splitting_moment_N3(0; log_power = 2, nf = nf),
            logarithmic_references.lo_second,
            rtol = 2.0e-9,
            atol = 2.0e-9,
        )
        test_moment_channels(
            splitting_moment_N3(1; log_power = 1, nf = nf),
            logarithmic_references.nlo_first,
            rtol = 2.0e-8,
            atol = 2.0e-7,
        )
    end
end

@testset "Momentum sum rule by loop coefficient" begin
    nf = 5
    tolerances = (1.0e-9, 1.0e-6, 1.0e-1)
    for loop_order in 0:2
        moment = splitting_moment(loop_order; mellin_N = 2, nf = nf)
        tolerance = tolerances[loop_order + 1]
        @test abs(moment.qq + moment.gq) < tolerance
        @test abs(moment.qg + moment.gg) < tolerance
    end
end
