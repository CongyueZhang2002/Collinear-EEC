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

@inline function diagonal_convolution_dot(
    shifted_vec::AbstractVector{Float64},
    current::Float64,
    kernel_vec::AbstractVector{Float64},
    d0_sub_weight_vec::AbstractVector{Float64},
    d1_sub_weight_vec::AbstractVector{Float64},
    d0_at_1::Float64,
    d1_at_1::Float64,
)

    if axes(shifted_vec, 1) != axes(kernel_vec, 1) ||
       axes(shifted_vec, 1) != axes(d0_sub_weight_vec, 1) ||
       axes(shifted_vec, 1) != axes(d1_sub_weight_vec, 1)
        throw(DimensionMismatch("Diagonal convolution vectors must have matching axes."))
    end

    result = 0.0
    @inbounds @simd for i in eachindex(
        shifted_vec,
        kernel_vec,
        d0_sub_weight_vec,
        d1_sub_weight_vec,
    )
        subtraction = (
            d0_at_1 * d0_sub_weight_vec[i] +
            d1_at_1 * d1_sub_weight_vec[i]
        )
        result += shifted_vec[i] * kernel_vec[i] - current * subtraction
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

    qq_closure_subtraction = (
        kernels.qq_d0_at_1_vec[time_index] *
        kernels.d0_closure_factor_vec[node_index] +
        kernels.qq_d1_at_1_vec[time_index] *
        kernels.d1_closure_factor_vec[node_index]
    )
    gg_closure_subtraction = (
        kernels.gg_d0_at_1_vec[time_index] *
        kernels.d0_closure_factor_vec[node_index] +
        kernels.gg_d1_at_1_vec[time_index] *
        kernels.d1_closure_factor_vec[node_index]
    )
    rhs_q = jq[node_index] * (
        kernels.qq_delta_at_1_vec[time_index] - qq_closure_subtraction
    )
    rhs_g = jg[node_index] * (
        kernels.gg_delta_at_1_vec[time_index] - gg_closure_subtraction
    )

    n_shift = node_index - 1
    if n_shift == 0
        return rhs_q, rhs_g
    end

    jq_shifted_vec = shifted_values(jq, node_index)
    jg_shifted_vec = shifted_values(jg, node_index)

    rhs_q += diagonal_convolution_dot(
        jq_shifted_vec,
        jq[node_index],
        @view(kernels.qq[time_index, 1:n_shift]),
        @view(kernels.d0_sub_weight_vec[1:n_shift]),
        @view(kernels.d1_sub_weight_vec[1:n_shift]),
        kernels.qq_d0_at_1_vec[time_index],
        kernels.qq_d1_at_1_vec[time_index],
    )
    rhs_q += convolution_dot(jg_shifted_vec, @view kernels.gq[time_index, 1:n_shift])

    rhs_g += convolution_dot(jq_shifted_vec, @view kernels.qg[time_index, 1:n_shift])
    rhs_g += diagonal_convolution_dot(
        jg_shifted_vec,
        jg[node_index],
        @view(kernels.gg[time_index, 1:n_shift]),
        @view(kernels.d0_sub_weight_vec[1:n_shift]),
        @view(kernels.d1_sub_weight_vec[1:n_shift]),
        kernels.gg_d0_at_1_vec[time_index],
        kernels.gg_d1_at_1_vec[time_index],
    )

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
