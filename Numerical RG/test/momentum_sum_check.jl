using Printf
using QuadGK

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

function splitting_moment_N2(; as::Float64, order::Int64, nf::Int64)
    function components(y)
        return timelike_splitting_convolution_func(
            y = y, as = as, order = order, nf = nf,
        )
    end

    sample = components(0.5)

    qq_integrand(y) = begin
        value = components(y)
        (
            y * value.PqqReg +
            (y * value.PqqD0 - sample.PqqD0_at_1) / (1.0 - y) +
            (y * value.PqqD1 - sample.PqqD1_at_1) *
            log1p(-y) / (1.0 - y)
        )
    end
    gg_integrand(y) = begin
        value = components(y)
        (
            y * value.PggReg +
            (y * value.PggD0 - sample.PggD0_at_1) / (1.0 - y) +
            (y * value.PggD1 - sample.PggD1_at_1) *
            log1p(-y) / (1.0 - y)
        )
    end

    qq = sample.PqqDelta_at_1 + quadgk(qq_integrand, 0.0, 1.0, rtol = 1e-9)[1]
    qg = quadgk(y -> y * components(y).Pqg, 0.0, 1.0, rtol = 1e-9)[1]
    gq = quadgk(y -> y * components(y).Pgq, 0.0, 1.0, rtol = 1e-9)[1]
    gg = sample.PggDelta_at_1 + quadgk(gg_integrand, 0.0, 1.0, rtol = 1e-9)[1]

    return (qq = qq, qg = qg, gq = gq, gg = gg)
end

const MOMENT_NF_SCHEME = :VFNS
const MOMENT_MU = 20.0
for order in 0:2
    moments = splitting_moment_N2(
        as = 0.01,
        order = order,
        nf = nf_func(MOMENT_MU; scheme = MOMENT_NF_SCHEME),
    )
    q_column_sum = moments.qq + moments.gq
    g_column_sum = moments.qg + moments.gg
    @printf(
        "MOMENT_SUM order=%d qq=%.12e qg=%.12e gq=%.12e gg=%.12e q_sum=%.12e g_sum=%.12e\n",
        order,
        moments.qq,
        moments.qg,
        moments.gq,
        moments.gg,
        q_column_sum,
        g_column_sum,
    )
end
