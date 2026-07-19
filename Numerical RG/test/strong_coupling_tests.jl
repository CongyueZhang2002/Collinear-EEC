using Test

if !isdefined(@__MODULE__, :alpha_s_grid)
    include(joinpath(@__DIR__, "..", "numerical_rg.jl"))
end

const LEGACY_BETA_COEFFICIENTS = (
    7.66667,
    38.6665,
    180.908,
    4826.18,
)

function legacy_beta_alpha_s(alpha_s::Float64; order::Int64 = 4)
    a_s = alpha_s / (4.0 * pi)
    series = 0.0
    a_s_power = a_s
    for coefficient in LEGACY_BETA_COEFFICIENTS[1:order]
        series += coefficient * a_s_power
        a_s_power *= a_s
    end
    return -2.0 * alpha_s * series
end

function legacy_alpha_s_numerical(
    mu::Float64;
    mu_ref::Float64 = 91.2,
    alpha_s_ref::Float64 = 0.118,
    order::Int64 = 4,
)
    rhs(alpha_s, _, _) = legacy_beta_alpha_s(alpha_s; order = order)
    problem = ODEProblem(rhs, alpha_s_ref, (log(mu_ref), log(mu)))
    solution = solve(problem, Tsit5(); reltol = 1.0e-11, abstol = 1.0e-13)
    return solution[end]
end

# This reproduces the old repository's public core/alpha_s.jl formula.
function legacy_alpha_s_resummed(
    mu::Float64;
    mu_ref::Float64 = 91.2,
    alpha_s_ref::Float64 = 0.118,
)
    beta0, beta1, beta2, beta3 = LEGACY_BETA_COEFFICIENTS
    ell = 1.0 + beta0 * alpha_s_ref * log(mu^2 / mu_ref^2) / (4.0 * pi)

    order2 = -alpha_s_ref / (4.0 * pi * ell) * beta1 / beta0 * log(ell)
    order3 = (alpha_s_ref / (4.0 * pi * ell))^2 * (
        (beta1 / beta0)^2 * (log(ell)^2 - log(ell) + ell - 1.0) -
        beta2 / beta0 * (ell - 1.0)
    )
    order4 = (alpha_s_ref / (4.0 * pi * ell))^3 * (
        (beta1 / beta0)^3 * (
            -log(ell)^3 + 2.5 * log(ell)^2 -
            2.0 * (ell - 1.0) * log(ell) - 0.5 * (ell - 1.0)^2
        ) +
        beta1 * beta2 / beta0^2 * (
            (ell - 1.0) * ell + (2.0 * ell - 3.0) * log(ell)
        ) +
        beta3 / beta0 * (1.0 - ell^2) / 2.0
    )
    return alpha_s_ref / ell * (1.0 + order2 + order3 + order4)
end

function one_loop_vfns_reference(mu::Float64)
    points = Float64[MZ]
    if mu < MZ
        for mass in (mb, mc)
            mu < mass < MZ && push!(points, mass)
        end
        push!(points, mu)
    else
        mu > mt > MZ && push!(points, mt)
        push!(points, mu)
    end

    inverse_alpha_s = 1.0 / DEFAULT_ALPHA_S_MZ
    for index in 1:(length(points) - 1)
        mu_left = points[index]
        mu_right = points[index + 1]
        mu_midpoint = sqrt(mu_left * mu_right)
        nf = nf_func(mu_midpoint; scheme = :VFNS)
        inverse_alpha_s += β0_func(nf) * log(mu_right / mu_left) / (2.0 * pi)
    end
    return 1.0 / inverse_alpha_s
end

@testset "Flavor schemes" begin
    @test NF_SCHEMES == (:VFNS, :nf5)
    @test nf_func(1.0; scheme = :VFNS) == 3
    @test nf_func(prevfloat(mc); scheme = :VFNS) == 3
    @test nf_func(mc; scheme = :VFNS) == 4
    @test nf_func(prevfloat(mb); scheme = :VFNS) == 4
    @test nf_func(mb; scheme = :VFNS) == 5
    @test nf_func(prevfloat(mt); scheme = :VFNS) == 5
    @test nf_func(mt; scheme = :VFNS) == 6
    @test nf_func(10.0) == 5

    for mu in (0.5, mc, mb, mt, 1.0e4)
        @test nf_func(mu; scheme = :nf5) == 5
    end

    for scheme in NF_SCHEMES, mu in (2.0, 10.0, 300.0)
        nf = nf_func(mu; scheme = scheme)
        @test β_func(
            αs = 0.2,
            order = 4,
            mu = mu,
            nf_scheme = scheme,
        ) == beta_alpha_s(0.2; order = 4, nf = nf)
    end

    @test_throws ArgumentError nf_func(10.0; scheme = :fixed4)
    @test_throws DomainError nf_func(0.0; scheme = :VFNS)
    @test_throws DomainError nf_func(Inf; scheme = :VFNS)
end

@testset "Numerical strong coupling" begin
    alpha_s = 0.2
    beta3_references = (
        1.209037813024e4,
        8.035186419171e3,
        4.826156328096e3,
    )
    beta4_references = (
        1.303779068020e5,
        5.831055395044e4,
        1.547061222594e4,
    )
    for (index, nf) in enumerate(3:5)
        @test β3_func(nf) ≈ beta3_references[index] rtol = 2.0e-10
        @test β4_func(nf) ≈ beta4_references[index] rtol = 2.0e-9
    end

    for nf in 3:6
        expected_five_loop_increment = (
            -2.0 * alpha_s * (alpha_s / (4.0 * pi))^5 * β4_func(nf)
        )
        @test (
            beta_alpha_s(alpha_s; order = 5, nf = nf) -
            beta_alpha_s(alpha_s; order = 4, nf = nf)
        ) ≈ expected_five_loop_increment
    end

    scales = [1.1, 2.0, mb, 10.0, MZ, mt, 300.0]
    for scheme in NF_SCHEMES
        at_reference = alpha_s_func_numerical(
            mu_f = MZ,
            order = 4,
            nf_scheme = scheme,
        )
        @test at_reference == DEFAULT_ALPHA_S_MZ

        values = alpha_s_grid(
            log.(scales.^2);
            order = 4,
            nf_scheme = scheme,
        )
        @test all(>(0.0), values)
        @test all(<(0.0), diff(values))
        @test values[5] == DEFAULT_ALPHA_S_MZ
        @test αs_func(
            10.0;
            order = 4,
            nf_scheme = scheme,
        ) ≈ values[4] rtol = 2.0e-9
        @test isfinite(Γ_func(mu = 10.0, order = 2, nf_scheme = scheme))
    end

    for scale in (1.1, 2.0, 10.0, 300.0)
        observed = alpha_s_func_numerical(
            mu_f = scale,
            order = 1,
            nf_scheme = :VFNS,
        )
        @test observed ≈ one_loop_vfns_reference(scale) rtol = 2.0e-8
    end

    fixed_scale = 10.0
    fixed_one_loop = alpha_s_func_numerical(
        mu_f = fixed_scale,
        order = 1,
        nf_scheme = :nf5,
    )
    expected_fixed_one_loop = DEFAULT_ALPHA_S_MZ / (
        1.0 + β0_func(5) * DEFAULT_ALPHA_S_MZ *
        log(fixed_scale / MZ) / (2.0 * pi)
    )
    @test fixed_one_loop ≈ expected_fixed_one_loop rtol = 2.0e-10

    for threshold in (mc, mb, mt)
        below = αs_func(
            threshold * (1.0 - 1.0e-8);
            order = 4,
            nf_scheme = :VFNS,
        )
        above = αs_func(
            threshold * (1.0 + 1.0e-8);
            order = 4,
            nf_scheme = :VFNS,
        )
        @test below ≈ above rtol = 1.0e-7
    end

    @test_throws ArgumentError alpha_s_grid(
        log.([10.0, 20.0].^2);
        order = 4,
        nf_scheme = :invalid,
    )
end

@testset "Strong coupling versus legacy implementation" begin
    scales = Float64[2.0, 4.18, 5.0, 10.0, 20.0, 50.0, 91.2, 100.0, 160.0, 500.0]
    new_values = [
        alpha_s_func_numerical(
            mu_f = mu,
            mu_ref = 91.2,
            order = 4,
            nf_scheme = :nf5,
        ) for mu in scales
    ]
    old_numerical_values = legacy_alpha_s_numerical.(scales)
    old_resummed_values = legacy_alpha_s_resummed.(scales)

    numerical_relative_errors = abs.(
        (new_values .- old_numerical_values) ./ old_numerical_values
    )
    resummed_relative_errors = abs.(
        (new_values .- old_resummed_values) ./ old_resummed_values
    )

    @test maximum(numerical_relative_errors) < 3.0e-7
    @test maximum(resummed_relative_errors) < 1.7e-3
    @test new_values[7] == old_numerical_values[7] == old_resummed_values[7]
end
