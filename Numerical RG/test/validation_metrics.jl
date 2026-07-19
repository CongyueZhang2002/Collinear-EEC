using Printf
using Test

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))
include(joinpath(@__DIR__, "strong_coupling_tests.jl"))
include(joinpath(@__DIR__, "manufactured_nonlocal_tests.jl"))
include(joinpath(@__DIR__, "physical_convolution_tests.jl"))
include(joinpath(@__DIR__, "paper_moment_tests.jl"))

const METRICS_NF_SCHEME = :VFNS
const METRICS_MU = 20.0

function print_fixed_nf5_coupling_differences()
    scales = Float64[
        2.0, 4.18, 5.0, 10.0, 20.0, 50.0, 91.2, 100.0, 160.0, 500.0,
    ]
    numerical_values = [
        alpha_s_func_numerical(
            mu_f = mu,
            mu_ref = 91.2,
            order = 4,
            nf_scheme = :nf5,
        ) for mu in scales
    ]
    rounded_ode_values = legacy_alpha_s_numerical.(scales)
    truncated_values = legacy_alpha_s_resummed.(scales)

    rounded_ode_difference = maximum(abs.(
        numerical_values ./ rounded_ode_values .- 1.0
    ))
    truncated_difference = maximum(abs.(
        numerical_values ./ truncated_values .- 1.0
    ))
    @printf(
        "ALPHA_NF5 rounded_ode_max_relative_difference=%.12e truncated_max_relative_difference=%.12e\n",
        rounded_ode_difference,
        truncated_difference,
    )
end

print_fixed_nf5_coupling_differences()

function maximum_relative_moment_difference(actual, reference)
    return maximum(
        abs(getproperty(actual, channel) / getproperty(reference, channel) - 1.0)
        for channel in propertynames(reference)
    )
end

function normalized_column_residual(first_entry, second_entry)
    scale = abs(first_entry) + abs(second_entry)
    return scale == 0.0 ? 0.0 : abs(first_entry + second_entry) / scale
end

for loop_order in 0:2
    maximum_relative_difference = 0.0
    for nf in 3:6
        reference = paper_N3_references(nf)[loop_order + 1]
        actual = splitting_moment_N3(loop_order; nf = nf)
        maximum_relative_difference = max(
            maximum_relative_difference,
            maximum_relative_moment_difference(actual, reference),
        )
    end
    @printf(
        "PUBLISHED_N3 loop=%d maximum_relative_difference=%.12e\n",
        loop_order,
        maximum_relative_difference,
    )
end

function maximum_logarithmic_moment_difference()
    maximum_relative_difference = 0.0
    for nf in 3:6
        reference = paper_log_moment_references(nf)
        comparisons = (
            (
                splitting_moment_N3(0; log_power = 1, nf = nf),
                reference.lo_first,
            ),
            (
                splitting_moment_N3(0; log_power = 2, nf = nf),
                reference.lo_second,
            ),
            (
                splitting_moment_N3(1; log_power = 1, nf = nf),
                reference.nlo_first,
            ),
        )
        for (actual, expected) in comparisons
            maximum_relative_difference = max(
                maximum_relative_difference,
                maximum_relative_moment_difference(actual, expected),
            )
        end
    end
    return maximum_relative_difference
end

@printf(
    "PUBLISHED_LOG maximum_relative_difference=%.12e\n",
    maximum_logarithmic_moment_difference(),
)

for method in (:euler, :rk2)
    coarse = manufactured_nonlocal_error(17, method)
    fine = manufactured_nonlocal_error(33, method)
    observed_order = log(coarse.relative_l2 / fine.relative_l2) /
                     log(coarse.maximum_delta_t / fine.maximum_delta_t)
    @printf(
        "NONLOCAL method=%s coarse_l2=%.12e fine_l2=%.12e order=%.6f\n",
        string(method),
        coarse.relative_l2,
        fine.relative_l2,
        observed_order,
    )
end

for order in 0:2
    medium = physical_convolution_error(129, Int64(order), METRICS_NF_SCHEME)
    fine = physical_convolution_error(257, Int64(order), METRICS_NF_SCHEME)
    q_order = log(medium.q_relative_error / fine.q_relative_error) /
              log(medium.delta_ell / fine.delta_ell)
    g_order = log(medium.g_relative_error / fine.g_relative_error) /
              log(medium.delta_ell / fine.delta_ell)
    @printf(
        "PHYSICAL order=%d q_error=%.12e g_error=%.12e q_order=%.6f g_order=%.6f\n",
        order,
        fine.q_relative_error,
        fine.g_relative_error,
        q_order,
        g_order,
    )
end

for loop_order in 0:2
    moment = splitting_moment(
        loop_order,
        mellin_N = 2,
        nf = nf_func(METRICS_MU; scheme = METRICS_NF_SCHEME),
    )
    @printf(
        "MOMENTUM loop=%d q_normalized_residual=%.12e g_normalized_residual=%.12e\n",
        loop_order,
        normalized_column_residual(moment.qq, moment.gq),
        normalized_column_residual(moment.qg, moment.gg),
    )
end

function print_cumulative_momentum_residuals()
    cumulative_moment = (qq = 0.0, qg = 0.0, gq = 0.0, gg = 0.0)
    for loop_order in 0:2
        coefficient = splitting_moment(
            loop_order,
            mellin_N = 2,
            nf = nf_func(METRICS_MU; scheme = METRICS_NF_SCHEME),
        )
        weight = 0.01^(loop_order + 1)
        cumulative_moment = (
            qq = cumulative_moment.qq + weight * coefficient.qq,
            qg = cumulative_moment.qg + weight * coefficient.qg,
            gq = cumulative_moment.gq + weight * coefficient.gq,
            gg = cumulative_moment.gg + weight * coefficient.gg,
        )
        @printf(
            "CUMULATIVE_MOMENTUM order=%d q_normalized_residual=%.12e g_normalized_residual=%.12e\n",
            loop_order,
            normalized_column_residual(cumulative_moment.qq, cumulative_moment.gq),
            normalized_column_residual(cumulative_moment.qg, cumulative_moment.gg),
        )
    end
end

print_cumulative_momentum_residuals()
