struct SplittingKernelTable
    qq_delta_at_1_vec::Vector{Float64}
    gg_delta_at_1_vec::Vector{Float64}
    qq::Matrix{Float64}
    qg::Matrix{Float64}
    gq::Matrix{Float64}
    gg::Matrix{Float64}
    qq_self_sub_vec::Vector{Float64}
    gg_self_sub_vec::Vector{Float64}
end

function set_alpha_s_runtime!(; alpha_s_Z_value::Float64, alpha_s_loops::Int64)
    if !isfinite(alpha_s_Z_value) || alpha_s_Z_value <= 0.0
        throw(ArgumentError("alpha_s_Z_value must be finite and positive."))
    end
    if !(1 <= alpha_s_loops <= 4)
        throw(ArgumentError("alpha_s_loops must be between 1 and 4."))
    end

    global αs_Z = alpha_s_Z_value
    global nloops_αs = alpha_s_loops
    return nothing
end

function build_splitting_kernel_table(;
    delta_ell::Float64,
    n_nodes::Int64,
    t_grid::AbstractVector{Float64},
    order::Int64,
    alpha_s_Z_value::Float64 = 0.118,
    alpha_s_loops::Int64 = 4,
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
    if !isfinite(alpha_s_max) || alpha_s_max <= 0.0
        throw(ArgumentError("alpha_s_max must be finite and positive."))
    end

    set_alpha_s_runtime!(alpha_s_Z_value = alpha_s_Z_value, alpha_s_loops = alpha_s_loops)

    n_shift = n_nodes - 1
    shift_vec = collect(1:n_shift)
    y_vec = exp.(-delta_ell .* shift_vec)
    shifted_weight_vec = delta_ell .* y_vec.^3
    subtraction_weight_vec = delta_ell .* y_vec
    one_minus_y_vec = 1 .- y_vec
    log_one_minus_y_vec = log.(one_minus_y_vec)
    y_ref = exp(-delta_ell)

    if !(0.0 < y_ref < 1.0)
        throw(DomainError(delta_ell, "delta_ell is outside the representable kernel-grid range."))
    end

    tail_log = log1p(-y_vec[end])
    d0_self_factor = sum(subtraction_weight_vec ./ one_minus_y_vec) - tail_log
    d1_self_factor = (
        sum(subtraction_weight_vec .* log_one_minus_y_vec ./ one_minus_y_vec) -
        0.5 * tail_log^2
    )

    qq_delta_at_1_vec = zeros(Float64, n_nodes)
    gg_delta_at_1_vec = zeros(Float64, n_nodes)
    qq = zeros(Float64, n_nodes, n_shift)
    qg = zeros(Float64, n_nodes, n_shift)
    gq = zeros(Float64, n_nodes, n_shift)
    gg = zeros(Float64, n_nodes, n_shift)
    qq_self_sub_vec = zeros(Float64, n_nodes)
    gg_self_sub_vec = zeros(Float64, n_nodes)
    a_s_vec = zeros(Float64, n_nodes)

    for n in 1:n_nodes
        mu = exp(0.5 * t_grid[n])
        alpha_s = try
            alpha_s_func(mu)
        catch error
            throw(DomainError(
                mu,
                "alpha_s_func failed at mu = $mu: $(sprint(showerror, error))",
            ))
        end
        if !isfinite(alpha_s) || alpha_s <= 0.0 || alpha_s > alpha_s_max
            throw(DomainError(
                alpha_s,
                "alpha_s(mu = $mu) must lie in (0, $alpha_s_max].",
            ))
        end
        a_s_vec[n] = alpha_s / (4 * pi)
    end

    Threads.@threads :static for n in 1:n_nodes
        a_s = a_s_vec[n]
        at_1_components = spliting_convolution_func(y = y_ref, as = a_s, order = order)

        qq_delta_at_1_vec[n] = at_1_components.PqqDelta_at_1
        gg_delta_at_1_vec[n] = at_1_components.PggDelta_at_1
        qq_self_sub_vec[n] = (
            at_1_components.PqqD0_at_1 * d0_self_factor +
            at_1_components.PqqD1_at_1 * d1_self_factor
        )
        gg_self_sub_vec[n] = (
            at_1_components.PggD0_at_1 * d0_self_factor +
            at_1_components.PggD1_at_1 * d1_self_factor
        )

        for k in 1:n_shift
            components = if k == 1
                at_1_components
            else
                spliting_convolution_func(y = y_vec[k], as = a_s, order = order)
            end

            qq[n, k] = shifted_weight_vec[k] * (
                components.PqqReg +
                components.PqqD0 / one_minus_y_vec[k] +
                components.PqqD1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )
            qg[n, k] = shifted_weight_vec[k] * components.Pqg

            gq[n, k] = shifted_weight_vec[k] * components.Pgq
            gg[n, k] = shifted_weight_vec[k] * (
                components.PggReg +
                components.PggD0 / one_minus_y_vec[k] +
                components.PggD1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )
        end
    end

    return SplittingKernelTable(
        qq_delta_at_1_vec,
        gg_delta_at_1_vec,
        qq,
        qg,
        gq,
        gg,
        qq_self_sub_vec,
        gg_self_sub_vec,
    )
end
