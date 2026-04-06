@inline function t_to_mu(t::Float64)
    return exp(0.5 * t)
end

@inline function shifted_values(values::AbstractVector{Float64}, node_index::Int64)
    if node_index <= 1
        return @view values[1:0]
    end
    return @view values[node_index-1:-1:1]
end

@inline function convolution_dot(values::AbstractVector{Float64}, kernel::AbstractVector{Float64})
    return dot(values, kernel)
end

function build_rhs(;
    jq::AbstractVector{Float64},
    jg::AbstractVector{Float64},
    active::AbstractVector{Bool},
    components::SplittingComponentGrid,
)

    n_nodes = length(jq)
    rhs_q = zeros(Float64, n_nodes)
    rhs_g = zeros(Float64, n_nodes)

    for i in 1:n_nodes
        if !active[i]
            continue
        end

        # The m = 0 contribution is carried by the explicit delta-at-1 terms.
        rhs_q[i] = jq[i] * components.PqqDelta_at_1
        rhs_g[i] = jg[i] * components.PggDelta_at_1

        n_shift = i - 1
        if n_shift == 0
            continue
        end

        q_shifted = shifted_values(jq, i)
        g_shifted = shifted_values(jg, i)

        rhs_q[i] += convolution_dot(q_shifted, @view components.q_from_q_kernel[1:n_shift])
        rhs_q[i] += convolution_dot(g_shifted, @view components.q_from_g_kernel[1:n_shift])
        rhs_q[i] -= jq[i] * components.q_self_sub[n_shift]

        rhs_g[i] += convolution_dot(q_shifted, @view components.g_from_q_kernel[1:n_shift])
        rhs_g[i] += convolution_dot(g_shifted, @view components.g_from_g_kernel[1:n_shift])
        rhs_g[i] -= jg[i] * components.g_self_sub[n_shift]
    end

    return (rhs_q = rhs_q, rhs_g = rhs_g)
end
