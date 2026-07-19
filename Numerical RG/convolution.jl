@inline function shifted_values(values_vec::AbstractVector{Float64}, node_index::Int64)
    if node_index <= 1
        return @view values_vec[1:0]
    end
    return @view values_vec[node_index-1:-1:1]
end

@inline function convolution_dot(values_vec::AbstractVector{Float64}, kernel_vec::AbstractVector{Float64})
    if axes(values_vec, 1) != axes(kernel_vec, 1)
        throw(DimensionMismatch("Convolution vectors must have matching axes."))
    end

    result = 0.0
    @inbounds @simd for i in eachindex(values_vec, kernel_vec)
        result += values_vec[i] * kernel_vec[i]
    end
    return result
end

const MAX_RHS_TASKS = 4
const MIN_PARALLEL_RHS_NODES = 128

@inline function _build_node_rhs(
    jq::AbstractVector{Float64},
    jg::AbstractVector{Float64},
    kernels::SplittingKernelTable,
    time_index::Int64,
    node_index::Int64,
)

    rhs_q = jq[node_index] * (
        kernels.qq_delta_at_1_vec[time_index] -
        kernels.qq_self_sub_vec[time_index]
    )
    rhs_g = jg[node_index] * (
        kernels.gg_delta_at_1_vec[time_index] -
        kernels.gg_self_sub_vec[time_index]
    )

    n_shift = node_index - 1
    if n_shift == 0
        return rhs_q, rhs_g
    end

    jq_shifted_vec = shifted_values(jq, node_index)
    jg_shifted_vec = shifted_values(jg, node_index)

    rhs_q += convolution_dot(jq_shifted_vec, @view kernels.qq[time_index, 1:n_shift])
    rhs_q += convolution_dot(jg_shifted_vec, @view kernels.gq[time_index, 1:n_shift])

    rhs_g += convolution_dot(jq_shifted_vec, @view kernels.qg[time_index, 1:n_shift])
    rhs_g += convolution_dot(jg_shifted_vec, @view kernels.gg[time_index, 1:n_shift])

    return rhs_q, rhs_g
end

function _fill_rhs_stride!(
    rhs_q::Vector{Float64},
    rhs_g::Vector{Float64},
    jq::AbstractVector{Float64},
    jg::AbstractVector{Float64},
    active::AbstractVector{Bool},
    kernels::SplittingKernelTable,
    time_index::Int64,
    first_node::Int64,
    last_node::Int64,
    stride::Int64,
)

    for i in first_node:stride:last_node
        if active[i]
            rhs_q[i], rhs_g[i] = _build_node_rhs(jq, jg, kernels, time_index, i)
        end
    end

    return nothing
end

function build_rhs(;
    jq::AbstractVector{Float64},
    jg::AbstractVector{Float64},
    active::AbstractVector{Bool},
    kernels::SplittingKernelTable,
    time_index::Int64,
)

    n_nodes = size(kernels.qq, 2) + 1
    expected_axis = Base.OneTo(n_nodes)

    if axes(jq, 1) != expected_axis ||
       axes(jg, 1) != expected_axis ||
       axes(active, 1) != expected_axis
        throw(DimensionMismatch("jq, jg, and active must match the kernel lattice."))
    end
    if !checkbounds(Bool, kernels.qq_delta_at_1_vec, time_index)
        throw(BoundsError(kernels.qq_delta_at_1_vec, time_index))
    end

    rhs_q = zeros(Float64, n_nodes)
    rhs_g = zeros(Float64, n_nodes)

    n_active = count(active)
    if n_active == 0
        return (rhs_q = rhs_q, rhs_g = rhs_g)
    end

    first_active = findfirst(active)::Int64
    last_active = findlast(active)::Int64
    n_tasks = min(MAX_RHS_TASKS, Threads.nthreads(), n_active)

    if n_tasks == 1 || n_active < MIN_PARALLEL_RHS_NODES
        _fill_rhs_stride!(
            rhs_q, rhs_g, jq, jg, active, kernels, time_index,
            first_active, last_active, 1,
        )
    else
        Threads.@sync for offset in 0:n_tasks-1
            first_node = first_active + offset
            Threads.@spawn _fill_rhs_stride!(
                rhs_q, rhs_g, jq, jg, active, kernels, time_index,
                first_node, last_active, n_tasks,
            )
        end
    end

    return (rhs_q = rhs_q, rhs_g = rhs_g)
end
