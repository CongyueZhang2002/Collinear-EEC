struct SplittingKernelTable
    qq_delta_at_1_vec::Vector{Float64}
    gg_delta_at_1_vec::Vector{Float64}
    qq::Matrix{Float64}
    qg::Matrix{Float64}
    gq::Matrix{Float64}
    gg::Matrix{Float64}
    qq_d0_at_1_vec::Vector{Float64}
    qq_d1_at_1_vec::Vector{Float64}
    gg_d0_at_1_vec::Vector{Float64}
    gg_d1_at_1_vec::Vector{Float64}
    d0_sub_weight_vec::Vector{Float64}
    d1_sub_weight_vec::Vector{Float64}
    d0_closure_factor_vec::Vector{Float64}
    d1_closure_factor_vec::Vector{Float64}
end

function build_splitting_kernel_table(;
    delta_ell::Float64,
    n_nodes::Int64,
    t_grid::AbstractVector{Float64},
    order::Int64,
    nf_scheme::Symbol = :VFNS,
    alpha_s_ref::Float64 = DEFAULT_ALPHA_S_MZ,
    alpha_s_mu_ref::Float64 = MZ,
    alpha_s_order::Int64 = 4,
    alpha_s_reltol::Float64 = 1.0e-10,
    alpha_s_abstol::Float64 = 1.0e-12,
    alpha_s_max::Float64 = 1.0,
)

    if n_nodes < 2
        throw(ArgumentError("n_nodes must be at least 2."))
    end
    if !isfinite(delta_ell) || delta_ell <= 0.0
        throw(ArgumentError("delta_ell must be finite and positive."))
    end
    if axes(t_grid, 1) != Base.OneTo(n_nodes)
        throw(DimensionMismatch("t_grid must have exactly n_nodes entries."))
    end
    if !all(isfinite, t_grid) || !all(>(0.0), diff(t_grid))
        throw(ArgumentError("t_grid must be finite and strictly increasing."))
    end
    if !(0 <= order <= 2)
        throw(ArgumentError("order must be 0, 1, or 2."))
    end
    _check_nf_scheme(nf_scheme)
    if !isfinite(alpha_s_max) || alpha_s_max <= 0.0
        throw(ArgumentError("alpha_s_max must be finite and positive."))
    end

    n_shift = n_nodes - 1
    shift_vec = collect(1:n_shift)
    y_vec = exp.(-delta_ell .* shift_vec)
    shifted_weight_vec = delta_ell .* y_vec.^3
    one_minus_y_vec = 1 .- y_vec
    log_one_minus_y_vec = log.(one_minus_y_vec)
    d0_sub_weight_vec = delta_ell .* y_vec ./ one_minus_y_vec
    d1_sub_weight_vec = d0_sub_weight_vec .* log_one_minus_y_vec
    y_ref = exp(-delta_ell)

    if !(0.0 < y_ref < 1.0)
        throw(DomainError(delta_ell, "delta_ell is outside the representable kernel-grid range."))
    end

    tail_log = log1p(-y_vec[end])
    d0_closure_factor_vec = zeros(Float64, n_nodes)
    d1_closure_factor_vec = zeros(Float64, n_nodes)
    d0_closure_factor_vec[end] = -tail_log
    d1_closure_factor_vec[end] = -0.5 * tail_log^2
    for i in (n_nodes - 1):-1:1
        d0_closure_factor_vec[i] = (
            d0_sub_weight_vec[i] + d0_closure_factor_vec[i + 1]
        )
        d1_closure_factor_vec[i] = (
            d1_sub_weight_vec[i] + d1_closure_factor_vec[i + 1]
        )
    end

    qq_delta_at_1_vec = zeros(Float64, n_nodes)
    gg_delta_at_1_vec = zeros(Float64, n_nodes)
    qq = zeros(Float64, n_nodes, n_shift)
    qg = zeros(Float64, n_nodes, n_shift)
    gq = zeros(Float64, n_nodes, n_shift)
    gg = zeros(Float64, n_nodes, n_shift)
    qq_d0_at_1_vec = zeros(Float64, n_nodes)
    qq_d1_at_1_vec = zeros(Float64, n_nodes)
    gg_d0_at_1_vec = zeros(Float64, n_nodes)
    gg_d1_at_1_vec = zeros(Float64, n_nodes)
    mu_vec = exp.(0.5 .* t_grid)
    nf_vec = Int64[nf_func(mu; scheme = nf_scheme) for mu in mu_vec]
    alpha_s_vec = alpha_s_grid(
        t_grid;
        mu_ref = alpha_s_mu_ref,
        alpha_s_ref = alpha_s_ref,
        order = alpha_s_order,
        nf_scheme = nf_scheme,
        reltol = alpha_s_reltol,
        abstol = alpha_s_abstol,
    )
    a_s_vec = similar(alpha_s_vec)
    for n in eachindex(alpha_s_vec)
        alpha_s = alpha_s_vec[n]
        if !isfinite(alpha_s) || alpha_s <= 0.0 || alpha_s > alpha_s_max
            throw(DomainError(
                alpha_s,
                "alpha_s at t_grid[$n] must lie in (0, $alpha_s_max].",
            ))
        end
        a_s_vec[n] = alpha_s / (4.0 * pi)
    end

    n_coefficients = order + 1
    n_flavors = 4
    qq_delta_coefficient = zeros(Float64, n_flavors, n_coefficients)
    gg_delta_coefficient = zeros(Float64, n_flavors, n_coefficients)
    qq_d0_at_1_coefficient = zeros(Float64, n_flavors, n_coefficients)
    qq_d1_at_1_coefficient = zeros(Float64, n_flavors, n_coefficients)
    gg_d0_at_1_coefficient = zeros(Float64, n_flavors, n_coefficients)
    gg_d1_at_1_coefficient = zeros(Float64, n_flavors, n_coefficients)
    qq_coefficient = zeros(Float64, n_flavors, n_coefficients, n_shift)
    qg_coefficient = zeros(Float64, n_flavors, n_coefficients, n_shift)
    gq_coefficient = zeros(Float64, n_flavors, n_coefficients, n_shift)
    gg_coefficient = zeros(Float64, n_flavors, n_coefficients, n_shift)

    flavor_is_cached = falses(n_flavors)
    for nf in nf_vec
        flavor_index = nf - 2
        flavor_is_cached[flavor_index] && continue
        flavor_is_cached[flavor_index] = true

        for coefficient_index in 1:n_coefficients
            loop_order = coefficient_index - 1
            endpoint = timelike_splitting_coefficient_func(
                y = y_ref,
                loop_order = loop_order,
                nf = nf,
            )
            qq_delta_coefficient[flavor_index, coefficient_index] =
                endpoint.PqqDelta_at_1
            gg_delta_coefficient[flavor_index, coefficient_index] =
                endpoint.PggDelta_at_1
            qq_d0_at_1_coefficient[flavor_index, coefficient_index] =
                endpoint.PqqD0_at_1
            qq_d1_at_1_coefficient[flavor_index, coefficient_index] =
                endpoint.PqqD1_at_1
            gg_d0_at_1_coefficient[flavor_index, coefficient_index] =
                endpoint.PggD0_at_1
            gg_d1_at_1_coefficient[flavor_index, coefficient_index] =
                endpoint.PggD1_at_1
        end

        Threads.@threads :static for k in 1:n_shift
            for coefficient_index in 1:n_coefficients
                components = timelike_splitting_coefficient_func(
                    y = y_vec[k],
                    loop_order = coefficient_index - 1,
                    nf = nf,
                )
                qq_coefficient[flavor_index, coefficient_index, k] =
                    shifted_weight_vec[k] * (
                        components.PqqReg +
                        components.PqqD0 / one_minus_y_vec[k] +
                        components.PqqD1 * log_one_minus_y_vec[k] /
                        one_minus_y_vec[k]
                    )
                qg_coefficient[flavor_index, coefficient_index, k] =
                    shifted_weight_vec[k] * components.Pqg
                gq_coefficient[flavor_index, coefficient_index, k] =
                    shifted_weight_vec[k] * components.Pgq
                gg_coefficient[flavor_index, coefficient_index, k] =
                    shifted_weight_vec[k] * (
                        components.PggReg +
                        components.PggD0 / one_minus_y_vec[k] +
                        components.PggD1 * log_one_minus_y_vec[k] /
                        one_minus_y_vec[k]
                    )
            end
        end
    end

    Threads.@threads :static for n in 1:n_nodes
        flavor_index = nf_vec[n] - 2
        a_s_power = a_s_vec[n]
        for coefficient_index in 1:n_coefficients
            qq_delta_at_1_vec[n] += (
                a_s_power * qq_delta_coefficient[flavor_index, coefficient_index]
            )
            gg_delta_at_1_vec[n] += (
                a_s_power * gg_delta_coefficient[flavor_index, coefficient_index]
            )
            qq_d0_at_1_vec[n] += (
                a_s_power *
                qq_d0_at_1_coefficient[flavor_index, coefficient_index]
            )
            qq_d1_at_1_vec[n] += (
                a_s_power *
                qq_d1_at_1_coefficient[flavor_index, coefficient_index]
            )
            gg_d0_at_1_vec[n] += (
                a_s_power *
                gg_d0_at_1_coefficient[flavor_index, coefficient_index]
            )
            gg_d1_at_1_vec[n] += (
                a_s_power *
                gg_d1_at_1_coefficient[flavor_index, coefficient_index]
            )

            @views qq[n, :] .+= (
                a_s_power .* qq_coefficient[flavor_index, coefficient_index, :]
            )
            @views qg[n, :] .+= (
                a_s_power .* qg_coefficient[flavor_index, coefficient_index, :]
            )
            @views gq[n, :] .+= (
                a_s_power .* gq_coefficient[flavor_index, coefficient_index, :]
            )
            @views gg[n, :] .+= (
                a_s_power .* gg_coefficient[flavor_index, coefficient_index, :]
            )

            a_s_power *= a_s_vec[n]
        end
    end

    return SplittingKernelTable(
        qq_delta_at_1_vec,
        gg_delta_at_1_vec,
        qq,
        qg,
        gq,
        gg,
        qq_d0_at_1_vec,
        qq_d1_at_1_vec,
        gg_d0_at_1_vec,
        gg_d1_at_1_vec,
        d0_sub_weight_vec,
        d1_sub_weight_vec,
        d0_closure_factor_vec,
        d1_closure_factor_vec,
    )
end
