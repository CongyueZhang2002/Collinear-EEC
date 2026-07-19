using Printf

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

const EXPERIMENT_MAX_RHS_TASKS = 4
const EXPERIMENT_MIN_PARALLEL_NODES = 128
const EXPERIMENT_NF_SCHEME = :VFNS

function convert_kernel_table(::Type{T}, table::SplittingKernelTable) where {T<:AbstractFloat}
    names = fieldnames(SplittingKernelTable)
    converted = map(name -> T.(getfield(table, name)), names)
    return NamedTuple{names}(converted)
end

@inline function experiment_shifted_values(values_vec, node_index::Int)
    return @view values_vec[node_index-1:-1:1]
end

@inline function experiment_dot(values_vec, kernel_vec)
    result = zero(eltype(values_vec))
    @inbounds @simd for i in eachindex(values_vec, kernel_vec)
        result += values_vec[i] * kernel_vec[i]
    end
    return result
end

@inline function experiment_diagonal_dot(
    shifted_vec,
    current,
    kernel_vec,
    d0_sub_weight_vec,
    d1_sub_weight_vec,
    d0_at_1,
    d1_at_1,
)
    result = zero(eltype(shifted_vec))
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

@inline function experiment_node_rhs(jq, jg, kernels, time_index::Int, node_index::Int)
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
    n_shift == 0 && return rhs_q, rhs_g

    jq_shifted_vec = experiment_shifted_values(jq, node_index)
    jg_shifted_vec = experiment_shifted_values(jg, node_index)
    rhs_q += experiment_diagonal_dot(
        jq_shifted_vec,
        jq[node_index],
        @view(kernels.qq[time_index, 1:n_shift]),
        @view(kernels.d0_sub_weight_vec[1:n_shift]),
        @view(kernels.d1_sub_weight_vec[1:n_shift]),
        kernels.qq_d0_at_1_vec[time_index],
        kernels.qq_d1_at_1_vec[time_index],
    )
    rhs_q += experiment_dot(
        jg_shifted_vec,
        @view(kernels.gq[time_index, 1:n_shift]),
    )
    rhs_g += experiment_dot(
        jq_shifted_vec,
        @view(kernels.qg[time_index, 1:n_shift]),
    )
    rhs_g += experiment_diagonal_dot(
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

function experiment_fill_rhs_stride!(
    rhs_q,
    rhs_g,
    jq,
    jg,
    active,
    kernels,
    time_index::Int,
    first_node::Int,
    last_node::Int,
    stride::Int,
)
    for i in first_node:stride:last_node
        if active[i]
            rhs_q[i], rhs_g[i] = experiment_node_rhs(
                jq,
                jg,
                kernels,
                time_index,
                i,
            )
        end
    end
    return nothing
end

function experiment_build_rhs(jq, jg, active, kernels, time_index::Int)
    n_nodes = length(jq)
    T = eltype(jq)
    rhs_q = zeros(T, n_nodes)
    rhs_g = zeros(T, n_nodes)
    n_active = count(active)
    n_active == 0 && return (rhs_q = rhs_q, rhs_g = rhs_g)

    first_active = findfirst(active)::Int
    last_active = findlast(active)::Int
    n_tasks = min(EXPERIMENT_MAX_RHS_TASKS, Threads.nthreads(), n_active)

    if n_tasks == 1 || n_active < EXPERIMENT_MIN_PARALLEL_NODES
        experiment_fill_rhs_stride!(
            rhs_q,
            rhs_g,
            jq,
            jg,
            active,
            kernels,
            time_index,
            first_active,
            last_active,
            1,
        )
    else
        Threads.@sync for offset in 0:(n_tasks - 1)
            first_node = first_active + offset
            Threads.@spawn experiment_fill_rhs_stride!(
                rhs_q,
                rhs_g,
                jq,
                jg,
                active,
                kernels,
                time_index,
                first_node,
                last_active,
                n_tasks,
            )
        end
    end
    return (rhs_q = rhs_q, rhs_g = rhs_g)
end

function experiment_fill_lower!(
    jq_history,
    jg_history,
    time_grid,
    kernels,
    midpoint_kernels,
)
    n_nodes = size(jq_history, 1)
    T = eltype(jq_history)
    half = T(0.5)
    jq_state = zeros(T, n_nodes)
    jg_state = zeros(T, n_nodes)
    jq_mid_state = zeros(T, n_nodes)
    jg_mid_state = zeros(T, n_nodes)
    active = falses(n_nodes)

    for n in n_nodes:-1:2
        dt = T(time_grid[n] - time_grid[n - 1])
        active[n] = true
        @views jq_state .= jq_history[:, n]
        @views jg_state .= jg_history[:, n]
        k1 = experiment_build_rhs(jq_state, jg_state, active, kernels, n)

        for i in 1:(n - 1)
            jq_mid_state[i] = half * (
                jq_history[i, n] + jq_history[i, n - 1]
            )
            jg_mid_state[i] = half * (
                jg_history[i, n] + jg_history[i, n - 1]
            )
        end
        for i in n:n_nodes
            jq_mid_state[i] = jq_state[i] - half * dt * k1.rhs_q[i]
            jg_mid_state[i] = jg_state[i] - half * dt * k1.rhs_g[i]
        end

        k2 = experiment_build_rhs(
            jq_mid_state,
            jg_mid_state,
            active,
            midpoint_kernels,
            n - 1,
        )
        for i in n:n_nodes
            jq_history[i, n - 1] = jq_state[i] - dt * k2.rhs_q[i]
            jg_history[i, n - 1] = jg_state[i] - dt * k2.rhs_g[i]
        end
    end
    return nothing
end

function experiment_evolve(
    ::Type{T},
    lattice::LatticeGrid,
    boundary_values,
    kernels,
    midpoint_kernels,
) where {T<:AbstractFloat}
    n_nodes = lattice.n_nodes
    half = T(0.5)
    jq_bc = T.(boundary_values.jq_bc)
    jg_bc = T.(boundary_values.jg_bc)
    jq = zeros(T, n_nodes)
    jg = zeros(T, n_nodes)
    active = falses(n_nodes)
    jq_history = zeros(T, n_nodes, n_nodes)
    jg_history = zeros(T, n_nodes, n_nodes)

    for n in 1:n_nodes
        jq[n] = jq_bc[n]
        jg[n] = jg_bc[n]
        active[n] = true
        jq_history[:, n] .= jq
        jg_history[:, n] .= jg
        n == n_nodes && break

        dt = T(lattice.t_grid[n + 1] - lattice.t_grid[n])
        k1 = experiment_build_rhs(jq, jg, active, kernels, n)
        jq_mid = copy(jq)
        jg_mid = copy(jg)
        for i in 1:n
            jq_mid[i] += half * dt * k1.rhs_q[i]
            jg_mid[i] += half * dt * k1.rhs_g[i]
        end
        k2 = experiment_build_rhs(
            jq_mid,
            jg_mid,
            active,
            midpoint_kernels,
            n,
        )
        for i in 1:n
            jq[i] += dt * k2.rhs_q[i]
            jg[i] += dt * k2.rhs_g[i]
        end
    end

    experiment_fill_lower!(
        jq_history,
        jg_history,
        lattice.t_grid,
        kernels,
        midpoint_kernels,
    )
    return (
        jq_history = jq_history,
        jg_history = jg_history,
        jq_bc = jq_bc,
        jg_bc = jg_bc,
    )
end

function build_rk2_tables(lattice::LatticeGrid, nf_scheme::Symbol)
    kernels = build_splitting_kernel_table(
        delta_ell = lattice.delta_ell,
        n_nodes = lattice.n_nodes,
        t_grid = lattice.t_grid,
        order = 2,
        nf_scheme = nf_scheme,
    )
    midpoint_grid = copy(lattice.t_grid)
    for n in 1:(lattice.n_nodes - 1)
        midpoint_grid[n] = 0.5 * (
            lattice.t_grid[n] + lattice.t_grid[n + 1]
        )
    end
    midpoint_grid[end] = lattice.t_grid[end]
    midpoint_kernels = build_splitting_kernel_table(
        delta_ell = lattice.delta_ell,
        n_nodes = lattice.n_nodes,
        t_grid = midpoint_grid,
        order = 2,
        nf_scheme = nf_scheme,
    )
    return kernels, midpoint_kernels
end

function best_timing(f; repeats::Int = 3)
    best = nothing
    for _ in 1:repeats
        GC.gc()
        measured = @timed f()
        if best === nothing || measured.time < best.time
            best = measured
        end
    end
    return best
end

function history_error(reference, candidate, lattice::LatticeGrid)
    q_scale = maximum(abs, reference.jq_history)
    g_scale = maximum(abs, reference.jg_history)
    q_cut = 1.0e-10 * q_scale
    g_cut = 1.0e-10 * g_scale
    q_error_sq = 0.0
    g_error_sq = 0.0
    q_reference_sq = 0.0
    g_reference_sq = 0.0
    upper_q_error_sq = 0.0
    upper_g_error_sq = 0.0
    upper_q_reference_sq = 0.0
    upper_g_reference_sq = 0.0
    lower_q_error_sq = 0.0
    lower_g_error_sq = 0.0
    lower_q_reference_sq = 0.0
    lower_g_reference_sq = 0.0
    max_absolute = 0.0
    max_significant_relative = 0.0
    column_relative = zeros(Float64, lattice.n_nodes)

    for n in 1:lattice.n_nodes
        column_error_sq = 0.0
        column_reference_sq = 0.0
        for i in 1:lattice.n_nodes
            q_ref = reference.jq_history[i, n]
            g_ref = reference.jg_history[i, n]
            q_error = Float64(candidate.jq_history[i, n]) - q_ref
            g_error = Float64(candidate.jg_history[i, n]) - g_ref
            q_error_sq += q_error^2
            g_error_sq += g_error^2
            q_reference_sq += q_ref^2
            g_reference_sq += g_ref^2
            if i <= n
                upper_q_error_sq += q_error^2
                upper_g_error_sq += g_error^2
                upper_q_reference_sq += q_ref^2
                upper_g_reference_sq += g_ref^2
            end
            if i >= n
                lower_q_error_sq += q_error^2
                lower_g_error_sq += g_error^2
                lower_q_reference_sq += q_ref^2
                lower_g_reference_sq += g_ref^2
            end
            column_error_sq += q_error^2 + g_error^2
            column_reference_sq += q_ref^2 + g_ref^2
            max_absolute = max(max_absolute, abs(q_error), abs(g_error))
            if abs(q_ref) > q_cut
                max_significant_relative = max(
                    max_significant_relative,
                    abs(q_error / q_ref),
                )
            end
            if abs(g_ref) > g_cut
                max_significant_relative = max(
                    max_significant_relative,
                    abs(g_error / g_ref),
                )
            end
        end
        column_relative[n] = sqrt(column_error_sq / column_reference_sq)
    end

    max_column_error, max_column_index = findmax(column_relative)
    diagonal_error = 0.0
    for i in 1:lattice.n_nodes
        diagonal_error = max(
            diagonal_error,
            abs(Float64(candidate.jq_history[i, i]) - reference.jq_history[i, i]),
            abs(Float64(candidate.jg_history[i, i]) - reference.jg_history[i, i]),
        )
    end

    return (
        q_relative_l2 = sqrt(q_error_sq / q_reference_sq),
        g_relative_l2 = sqrt(g_error_sq / g_reference_sq),
        upper_q_relative_l2 = sqrt(upper_q_error_sq / upper_q_reference_sq),
        upper_g_relative_l2 = sqrt(upper_g_error_sq / upper_g_reference_sq),
        lower_q_relative_l2 = sqrt(lower_q_error_sq / lower_q_reference_sq),
        lower_g_relative_l2 = sqrt(lower_g_error_sq / lower_g_reference_sq),
        max_absolute = max_absolute,
        max_significant_relative = max_significant_relative,
        diagonal_absolute = diagonal_error,
        first_column_relative = column_relative[1],
        middle_column_relative = column_relative[cld(lattice.n_nodes, 2)],
        last_column_relative = column_relative[end],
        max_column_relative = max_column_error,
        max_column_mu = exp(0.5 * lattice.t_grid[max_column_index]),
    )
end

function run_experiment(n_nodes::Int)
    bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)
    boundary_func(; b, bstar, mu_start) = (
        exp(-2.0 * b),
        0.5 * exp(-3.0 * b),
    )
    lattice = build_lattice_grid(
        n_nodes = n_nodes,
        b_min = 0.001,
        b_max = 30.0,
        bstar_func = bstar_func,
    )
    boundary_values = build_boundary_values(
        grid = lattice,
        boundary_func = boundary_func,
    )

    table_timing = best_timing(
        () -> build_rk2_tables(lattice, EXPERIMENT_NF_SCHEME),
    )
    kernels64, midpoint64 = table_timing.value
    conversion_timing = best_timing(
        () -> (
            convert_kernel_table(Float32, kernels64),
            convert_kernel_table(Float32, midpoint64),
        ),
    )
    kernels32, midpoint32 = conversion_timing.value

    evolution64 = best_timing(
        () -> experiment_evolve(
            Float64,
            lattice,
            boundary_values,
            kernels64,
            midpoint64,
        ),
    )
    evolution32 = best_timing(
        () -> experiment_evolve(
            Float32,
            lattice,
            boundary_values,
            kernels32,
            midpoint32,
        ),
    )
    kernel32_state64 = best_timing(
        () -> experiment_evolve(
            Float64,
            lattice,
            boundary_values,
            kernels32,
            midpoint32,
        ),
    )
    error = history_error(evolution64.value, evolution32.value, lattice)
    hybrid_error = history_error(
        evolution64.value,
        kernel32_state64.value,
        lattice,
    )

    table64_mib = Base.summarysize((kernels64, midpoint64)) / 1024.0^2
    table32_mib = Base.summarysize((kernels32, midpoint32)) / 1024.0^2
    history64_mib = Base.summarysize(evolution64.value) / 1024.0^2
    history32_mib = Base.summarysize(evolution32.value) / 1024.0^2

    middle = cld(n_nodes, 2)
    jq64 = boundary_values.jq_bc
    jg64 = boundary_values.jg_bc
    jq32 = Float32.(jq64)
    jg32 = Float32.(jg64)
    rhs64 = experiment_node_rhs(jq64, jg64, kernels64, middle, middle)
    rhs32 = experiment_node_rhs(jq32, jg32, kernels32, middle, middle)
    rhs_q_relative = abs(Float64(rhs32[1]) - rhs64[1]) / abs(rhs64[1])
    rhs_g_relative = abs(Float64(rhs32[2]) - rhs64[2]) / abs(rhs64[2])

    @printf("\nNODES %d\n", n_nodes)
    @printf(
        "TIME kernel64=%.6f conversion32=%.6f evolve64=%.6f evolve32=%.6f speedup=%.3f\n",
        table_timing.time,
        conversion_timing.time,
        evolution64.time,
        evolution32.time,
        evolution64.time / evolution32.time,
    )
    @printf(
        "HYBRID_TIME kernel32_state64=%.6f speedup=%.3f\n",
        kernel32_state64.time,
        evolution64.time / kernel32_state64.time,
    )
    @printf(
        "ALLOC evolve64_mib=%.3f evolve32_mib=%.3f ratio=%.3f\n",
        evolution64.bytes / 1024.0^2,
        evolution32.bytes / 1024.0^2,
        evolution32.bytes / evolution64.bytes,
    )
    @printf(
        "STORAGE table64_mib=%.3f table32_mib=%.3f history64_mib=%.3f history32_mib=%.3f\n",
        table64_mib,
        table32_mib,
        history64_mib,
        history32_mib,
    )
    @printf(
        "STORAGE_TOTAL float64_mib=%.3f float32_mib=%.3f hybrid_mib=%.3f\n",
        table64_mib + history64_mib,
        table32_mib + history32_mib,
        table32_mib + history64_mib,
    )
    @printf(
        "ERROR q_rel_l2=%.6e g_rel_l2=%.6e max_abs=%.6e max_significant_rel=%.6e diagonal_abs=%.6e\n",
        error.q_relative_l2,
        error.g_relative_l2,
        error.max_absolute,
        error.max_significant_relative,
        error.diagonal_absolute,
    )
    @printf(
        "ERROR_REGION upper_q=%.6e upper_g=%.6e lower_q=%.6e lower_g=%.6e finite32=%s\n",
        error.upper_q_relative_l2,
        error.upper_g_relative_l2,
        error.lower_q_relative_l2,
        error.lower_g_relative_l2,
        string(
            all(isfinite, evolution32.value.jq_history) &&
            all(isfinite, evolution32.value.jg_history),
        ),
    )
    @printf(
        "ERROR_BY_MU first=%.6e middle=%.6e last=%.6e max=%.6e at_mu=%.6g\n",
        error.first_column_relative,
        error.middle_column_relative,
        error.last_column_relative,
        error.max_column_relative,
        error.max_column_mu,
    )
    @printf(
        "HYBRID_ERROR q_rel_l2=%.6e g_rel_l2=%.6e upper_g=%.6e lower_g=%.6e max_significant_rel=%.6e\n",
        hybrid_error.q_relative_l2,
        hybrid_error.g_relative_l2,
        hybrid_error.upper_g_relative_l2,
        hybrid_error.lower_g_relative_l2,
        hybrid_error.max_significant_relative,
    )
    @printf(
        "RHS q_relative=%.6e g_relative=%.6e\n",
        rhs_q_relative,
        rhs_g_relative,
    )
    return nothing
end

println("Julia version: ", VERSION)
println("Julia threads: ", Threads.nthreads())
println("Experiment: Float64 grid/coupling/kernel evaluation; Float32 kernel storage and evolution")

# Compile both precision paths before collecting timings.
bstar_warmup(b) = b
boundary_warmup(; b, bstar, mu_start) = (exp(-b), 0.5 * exp(-b))
warmup_lattice = build_lattice_grid(
    n_nodes = 12,
    b_min = 0.1,
    b_max = 1.0,
    bstar_func = bstar_warmup,
)
warmup_boundary = build_boundary_values(
    grid = warmup_lattice,
    boundary_func = boundary_warmup,
)
warmup64 = build_rk2_tables(warmup_lattice, EXPERIMENT_NF_SCHEME)
warmup32 = map(table -> convert_kernel_table(Float32, table), warmup64)
warmup_experiment64 = experiment_evolve(
    Float64,
    warmup_lattice,
    warmup_boundary,
    warmup64...,
)
warmup_production = _evolve_stepwise_rg(
    lattice = warmup_lattice,
    boundary_values = warmup_boundary,
    kernel_table = warmup64[1],
    midpoint_kernel_table = warmup64[2],
    method = :rk2,
    store_history = true,
)
reference_difference = max(
    maximum(abs, warmup_experiment64.jq_history .- warmup_production.jq_history),
    maximum(abs, warmup_experiment64.jg_history .- warmup_production.jg_history),
)
reference_difference == 0.0 || error(
    "Experimental Float64 path differs from production by $reference_difference.",
)
experiment_evolve(Float32, warmup_lattice, warmup_boundary, warmup32...)
println("Reference check: experimental Float64 evolution exactly matches production")
warmup64 = nothing
warmup32 = nothing
GC.gc()

node_counts = isempty(ARGS) ? [500, 1000] : parse.(Int, ARGS)
for n_nodes in node_counts
    run_experiment(n_nodes)
end
