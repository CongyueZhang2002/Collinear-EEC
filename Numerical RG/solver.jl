function unit_boundary_value(; b::Float64, bstar::Float64, mu_start::Float64)
    return 1.0, 1.0
end

function activate_nodes!(
    jq::Vector{Float64},
    jg::Vector{Float64},
    active::AbstractVector{Bool},
    jq_bc::AbstractVector{Float64},
    jg_bc::AbstractVector{Float64},
    mu_i_grid::AbstractVector{Float64},
    t::Float64;
    activation_tol::Float64 = 1e-12,
)

    for i in eachindex(active)
        if !active[i] && t + activation_tol >= log(mu_i_grid[i]^2)
            jq[i] = jq_bc[i]
            jg[i] = jg_bc[i]
            active[i] = true
        end
    end

    return nothing
end

function fill_lower_triangle!(
    jq_history::Matrix{Float64},
    jg_history::Matrix{Float64},
    time_grid::Vector{Float64},
    kernels::SplittingKernelTable,
    midpoint_kernels::Union{Nothing,SplittingKernelTable},
    method::Symbol,
    n_nodes::Int64,
)

    jq_state = zeros(Float64, n_nodes)
    jg_state = zeros(Float64, n_nodes)
    jq_mid_state = zeros(Float64, n_nodes)
    jg_mid_state = zeros(Float64, n_nodes)
    lower_active = falses(n_nodes)

    for n in n_nodes:-1:2
        dt = time_grid[n] - time_grid[n - 1]
        lower_active[n] = true

        @views jq_state .= jq_history[:, n]
        @views jg_state .= jg_history[:, n]

        k1 = build_rhs(
            jq = jq_state,
            jg = jg_state,
            active = lower_active,
            kernels = kernels,
            time_index = n,
        )

        if method == :euler
            for i in n:n_nodes
                jq_history[i, n - 1] = jq_state[i] - dt * k1.rhs_q[i]
                jg_history[i, n - 1] = jg_state[i] - dt * k1.rhs_g[i]
            end
        else
            for i in 1:(n - 1)
                jq_mid_state[i] = 0.5 * (
                    jq_history[i, n] +
                    jq_history[i, n - 1]
                )
                jg_mid_state[i] = 0.5 * (
                    jg_history[i, n] +
                    jg_history[i, n - 1]
                )
            end
            for i in n:n_nodes
                jq_mid_state[i] = jq_state[i] - 0.5 * dt * k1.rhs_q[i]
                jg_mid_state[i] = jg_state[i] - 0.5 * dt * k1.rhs_g[i]
            end

            k2 = build_rhs(
                jq = jq_mid_state,
                jg = jg_mid_state,
                active = lower_active,
                kernels = midpoint_kernels,
                time_index = n - 1,
            )

            for i in n:n_nodes
                jq_history[i, n - 1] = jq_state[i] - dt * k2.rhs_q[i]
                jg_history[i, n - 1] = jg_state[i] - dt * k2.rhs_g[i]
            end
        end

        for i in n:n_nodes
            if !isfinite(jq_history[i, n - 1]) || !isfinite(jg_history[i, n - 1])
                throw(DomainError(
                    (jq_history[i, n - 1], jg_history[i, n - 1]),
                    "Backward evolution became non-finite at b node $i and time node $(n - 1).",
                ))
            end
        end
    end

    return nothing
end

function solve_stepwise_rg(;
    lattice::LatticeGrid,
    boundary_func::Function = unit_boundary_value,
    order::Int64,
    alpha_s_Z_value::Float64 = 0.118,
    alpha_s_loops::Int64 = 4,
    alpha_s_max::Float64 = 1.0,
    method::Symbol = :euler,
    store_history::Bool = true,
)

    if !(method in (:euler, :rk2))
        throw(ArgumentError("Unknown method: $method. Expected :euler or :rk2."))
    end

    boundary_values = build_boundary_values(grid = lattice, boundary_func = boundary_func)
    time_grid = lattice.t_grid
    kernel_table = build_splitting_kernel_table(
        delta_ell = lattice.delta_ell,
        n_nodes = lattice.n_nodes,
        t_grid = time_grid,
        order = order,
        alpha_s_Z_value = alpha_s_Z_value,
        alpha_s_loops = alpha_s_loops,
        alpha_s_max = alpha_s_max,
    )

    if method == :rk2
        midpoint_grid = copy(time_grid)

        for n in 1:(lattice.n_nodes - 1)
            midpoint_grid[n] = 0.5 * (time_grid[n] + time_grid[n + 1])
        end

        midpoint_grid[end] = time_grid[end]

        midpoint_kernel_table = build_splitting_kernel_table(
            delta_ell = lattice.delta_ell,
            n_nodes = lattice.n_nodes,
            t_grid = midpoint_grid,
            order = order,
            alpha_s_Z_value = alpha_s_Z_value,
            alpha_s_loops = alpha_s_loops,
            alpha_s_max = alpha_s_max,
        )
    elseif method == :euler
        midpoint_kernel_table = nothing
    end

    return _evolve_stepwise_rg(
        lattice = lattice,
        boundary_values = boundary_values,
        kernel_table = kernel_table,
        midpoint_kernel_table = midpoint_kernel_table,
        method = method,
        store_history = store_history,
    )
end

function _evolve_stepwise_rg(;
    lattice::LatticeGrid,
    boundary_values,
    kernel_table::SplittingKernelTable,
    midpoint_kernel_table::Union{Nothing,SplittingKernelTable},
    method::Symbol,
    store_history::Bool,
)

    if !(method in (:euler, :rk2))
        throw(ArgumentError("Unknown method: $method. Expected :euler or :rk2."))
    end
    if method == :rk2 && midpoint_kernel_table === nothing
        throw(ArgumentError("RK2 evolution requires a midpoint kernel table."))
    end

    time_grid = lattice.t_grid
    jq = zeros(Float64, lattice.n_nodes)
    jg = zeros(Float64, lattice.n_nodes)
    active = falses(lattice.n_nodes)

    if store_history
        jq_history = zeros(Float64, lattice.n_nodes, lattice.n_nodes)
        jg_history = zeros(Float64, lattice.n_nodes, lattice.n_nodes)
    else
        jq_history = zeros(Float64, 0, 0)
        jg_history = zeros(Float64, 0, 0)
    end

    for n in 1:lattice.n_nodes
        t_now = time_grid[n]

        activate_nodes!(
            jq,
            jg,
            active,
            boundary_values.jq_bc,
            boundary_values.jg_bc,
            lattice.mu_i_grid,
            t_now,
        )

        if store_history
            jq_history[:, n] .= jq
            jg_history[:, n] .= jg
        end

        if n == lattice.n_nodes
            break
        end

        dt = time_grid[n + 1] - t_now

        if method == :euler
            rhs = build_rhs(
                jq = jq,
                jg = jg,
                active = active,
                kernels = kernel_table,
                time_index = n,
            )

            for i in eachindex(active)
                if active[i]
                    jq[i] += dt * rhs.rhs_q[i]
                    jg[i] += dt * rhs.rhs_g[i]
                end
            end

        elseif method == :rk2
            k1 = build_rhs(
                jq = jq,
                jg = jg,
                active = active,
                kernels = kernel_table,
                time_index = n,
            )

            jq_mid = copy(jq)
            jg_mid = copy(jg)

            for i in eachindex(active)
                if active[i]
                    jq_mid[i] += 0.5 * dt * k1.rhs_q[i]
                    jg_mid[i] += 0.5 * dt * k1.rhs_g[i]
                end
            end

            k2 = build_rhs(
                jq = jq_mid,
                jg = jg_mid,
                active = active,
                kernels = midpoint_kernel_table,
                time_index = n,
            )

            for i in eachindex(active)
                if active[i]
                    jq[i] += dt * k2.rhs_q[i]
                    jg[i] += dt * k2.rhs_g[i]
                end
            end

        end
    end

    if store_history
        fill_lower_triangle!(
            jq_history,
            jg_history,
            time_grid,
            kernel_table,
            midpoint_kernel_table,
            method,
            lattice.n_nodes,
        )
    end

    return StepwiseRGSolution(
        jq,
        jg,
        active,
        boundary_values.jq_bc,
        boundary_values.jg_bc,
        jq_history,
        jg_history,
        lattice,
        collect(time_grid),
    )
end

"""
    solve_jet_rg(; n_nodes, b_min, b_max, bstar_func, boundary_func, order, ...)

Solve on the full rectangle bounded by the endpoint scales of
`mu_start(b) = b0 / bstar_func(b)`. `boundary_func(; b, bstar, mu_start)`
must return `(jq, jg)`. The returned `StepwiseRGSolution` is callable as
`solution(b, mu)`.

The convolution uses zero closure for `b > b_max`, so `b_max` must be large
enough that the physical jet function is negligible.
"""
function solve_jet_rg(;
    n_nodes::Int64,
    b_min::Float64,
    b_max::Float64,
    bstar_func::Function,
    boundary_func::Function = unit_boundary_value,
    order::Int64,
    alpha_s_Z_value::Float64 = 0.118,
    alpha_s_loops::Int64 = 4,
    alpha_s_max::Float64 = 1.0,
    method::Symbol = :rk2,
)

    lattice = build_lattice_grid(
        n_nodes = n_nodes,
        b_min = b_min,
        b_max = b_max,
        bstar_func = bstar_func,
    )

    return solve_stepwise_rg(
        lattice = lattice,
        boundary_func = boundary_func,
        order = order,
        alpha_s_Z_value = alpha_s_Z_value,
        alpha_s_loops = alpha_s_loops,
        alpha_s_max = alpha_s_max,
        method = method,
        store_history = true,
    )
end
