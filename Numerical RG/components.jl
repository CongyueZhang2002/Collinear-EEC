struct SplittingKernelTable
    qq_delta_at_1_vec::Vector{Float64}
    gg_delta_at_1_vec::Vector{Float64}
    qq::Matrix{Float64}
    qg::Matrix{Float64}
    gq::Matrix{Float64}
    gg::Matrix{Float64}
    qq_self_sub::Matrix{Float64}
    gg_self_sub::Matrix{Float64}
end

function set_alpha_s_runtime!(; alpha_s_Z_value::Float64, alpha_s_loops::Int64)
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
)

    set_alpha_s_runtime!(alpha_s_Z_value = alpha_s_Z_value, alpha_s_loops = alpha_s_loops)

    n_shift = n_nodes - 1
    shift_vec = collect(1:n_shift)
    y_vec = exp.(-delta_ell .* shift_vec)
    weight_vec = delta_ell .* y_vec.^3
    one_minus_y_vec = 1 .- y_vec
    log_one_minus_y_vec = log.(one_minus_y_vec)
    y_ref = exp(-delta_ell)

    qq_delta_at_1_vec = zeros(Float64, n_nodes)
    gg_delta_at_1_vec = zeros(Float64, n_nodes)
    qq = zeros(Float64, n_nodes, n_shift)
    qg = zeros(Float64, n_nodes, n_shift)
    gq = zeros(Float64, n_nodes, n_shift)
    gg = zeros(Float64, n_nodes, n_shift)
    qq_self_sub = zeros(Float64, n_nodes, n_shift)
    gg_self_sub = zeros(Float64, n_nodes, n_shift)

    for n in 1:n_nodes
        mu = exp(0.5 * t_grid[n])
        a_s = alpha_s_func(mu) / (4 * pi)
        at_1_components = spliting_convolution_func(y = y_ref, as = a_s, order = order)

        qq_delta_at_1_vec[n] = at_1_components.PqqDelta_at_1
        gg_delta_at_1_vec[n] = at_1_components.PggDelta_at_1

        q_self_running = 0.0
        g_self_running = 0.0

        for k in 1:n_shift
            components = spliting_convolution_func(y = y_vec[k], as = a_s, order = order)

            qq[n, k] = weight_vec[k] * (
                components.PqqReg +
                components.PqqD0 / one_minus_y_vec[k] +
                components.PqqD1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )
            qg[n, k] = weight_vec[k] * components.Pgq

            gq[n, k] = weight_vec[k] * components.Pqg
            gg[n, k] = weight_vec[k] * (
                components.PggReg +
                components.PggD0 / one_minus_y_vec[k] +
                components.PggD1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )

            q_self_running += weight_vec[k] * (
                at_1_components.PqqD0_at_1 / one_minus_y_vec[k] +
                at_1_components.PqqD1_at_1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )
            g_self_running += weight_vec[k] * (
                at_1_components.PggD0_at_1 / one_minus_y_vec[k] +
                at_1_components.PggD1_at_1 * log_one_minus_y_vec[k] / one_minus_y_vec[k]
            )

            qq_self_sub[n, k] = q_self_running
            gg_self_sub[n, k] = g_self_running
        end
    end

    return SplittingKernelTable(
        qq_delta_at_1_vec,
        gg_delta_at_1_vec,
        qq,
        qg,
        gq,
        gg,
        qq_self_sub,
        gg_self_sub,
    )
end
