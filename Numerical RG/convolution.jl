@inline function shifted_values(values_vec::AbstractVector{Float64}, node_index::Int64)
    if node_index <= 1
        return @view values_vec[1:0]
    end
    return @view values_vec[node_index-1:-1:1]
end

@inline function convolution_dot(values_vec::AbstractVector{Float64}, kernel_vec::AbstractVector{Float64})
    return dot(values_vec, kernel_vec)
end

function build_rhs(;
    jq::AbstractVector{Float64},
    jg::AbstractVector{Float64},
    active::AbstractVector{Bool},
    kernels::SplittingKernelTable,
    time_index::Int64,
)

    n_nodes = length(jq)
    rhs_q = zeros(Float64, n_nodes)
    rhs_g = zeros(Float64, n_nodes)

    for i in 1:n_nodes
        if !active[i]
            continue
        end

        # The m = 0 contribution is carried by the explicit delta-at-1 terms.
        rhs_q[i] = jq[i] * kernels.qq_delta_at_1_vec[time_index]
        rhs_g[i] = jg[i] * kernels.gg_delta_at_1_vec[time_index]

        n_shift = i - 1
        if n_shift == 0
            continue
        end

        jq_shifted_vec = shifted_values(jq, i)
        jg_shifted_vec = shifted_values(jg, i)

        rhs_q[i] += convolution_dot(jq_shifted_vec, @view kernels.qq[time_index, 1:n_shift])
        rhs_q[i] += convolution_dot(jg_shifted_vec, @view kernels.qg[time_index, 1:n_shift])
        rhs_q[i] -= jq[i] * kernels.qq_self_sub[time_index, n_shift]

        rhs_g[i] += convolution_dot(jq_shifted_vec, @view kernels.gq[time_index, 1:n_shift])
        rhs_g[i] += convolution_dot(jg_shifted_vec, @view kernels.gg[time_index, 1:n_shift])
        rhs_g[i] -= jg[i] * kernels.gg_self_sub[time_index, n_shift]
    end

    return (rhs_q = rhs_q, rhs_g = rhs_g)
end
