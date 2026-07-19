using LinearAlgebra

function manufactured_operator(n_nodes::Int64)
    operator = zeros(Float64, 2 * n_nodes, 2 * n_nodes)
    q_local = 0.12
    g_local = -0.08
    qq_shift = 0.07
    qg_shift = -0.05
    gq_shift = 0.09
    gg_shift = 0.04

    for i in 1:n_nodes
        operator[i, i] = q_local
        operator[n_nodes + i, n_nodes + i] = g_local
        if i > 1
            shifted_index = i - 1
            operator[i, shifted_index] = qq_shift
            operator[i, n_nodes + shifted_index] = gq_shift
            operator[n_nodes + i, shifted_index] = qg_shift
            operator[n_nodes + i, n_nodes + shifted_index] = gg_shift
        end
    end

    rates = (
        q_local = q_local,
        g_local = g_local,
        qq_shift = qq_shift,
        qg_shift = qg_shift,
        gq_shift = gq_shift,
        gg_shift = gg_shift,
    )
    return operator, rates
end

function manufactured_kernel_table(
    time_grid::Vector{Float64},
    rates;
    time_origin::Float64,
    time_slope::Float64,
)

    n_nodes = length(time_grid)
    n_shift = n_nodes - 1
    scale_vec = 1.0 .+ time_slope .* (time_grid .- time_origin)
    qq = zeros(Float64, n_nodes, n_shift)
    qg = zeros(Float64, n_nodes, n_shift)
    gq = zeros(Float64, n_nodes, n_shift)
    gg = zeros(Float64, n_nodes, n_shift)

    qq[:, 1] .= scale_vec .* rates.qq_shift
    qg[:, 1] .= scale_vec .* rates.qg_shift
    gq[:, 1] .= scale_vec .* rates.gq_shift
    gg[:, 1] .= scale_vec .* rates.gg_shift

    return SplittingKernelTable(
        scale_vec .* rates.q_local,
        scale_vec .* rates.g_local,
        qq,
        qg,
        gq,
        gg,
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
        zeros(Float64, n_shift),
        zeros(Float64, n_shift),
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
    )
end

function manufactured_nonlocal_error(n_nodes::Int64, method::Symbol)
    bstar_func(b) = b / sqrt(1.0 + (b / 0.8)^2)
    lattice = build_lattice_grid(
        n_nodes = n_nodes,
        b_min = 0.2,
        b_max = 2.0,
        bstar_func = bstar_func,
    )
    operator, rates = manufactured_operator(n_nodes)
    time_origin = lattice.t_grid[1]
    time_slope = 0.15
    initial_state = vcat(
        1.0 .+ 0.2 .* lattice.b_grid ./ lattice.b_max,
        0.8 .- 0.1 .* lattice.b_grid ./ lattice.b_max,
    )

    exact_q = zeros(Float64, n_nodes, n_nodes)
    exact_g = zeros(Float64, n_nodes, n_nodes)
    for n in 1:n_nodes
        delta_t = lattice.t_grid[n] - time_origin
        integrated_scale = delta_t + 0.5 * time_slope * delta_t^2
        state = exp(integrated_scale * operator) * initial_state
        exact_q[:, n] .= @view state[1:n_nodes]
        exact_g[:, n] .= @view state[(n_nodes + 1):(2 * n_nodes)]
    end

    boundary_values = (
        jq_bc = [exact_q[i, i] for i in 1:n_nodes],
        jg_bc = [exact_g[i, i] for i in 1:n_nodes],
    )
    kernel_table = manufactured_kernel_table(
        lattice.t_grid,
        rates,
        time_origin = time_origin,
        time_slope = time_slope,
    )

    midpoint_kernel_table = if method == :rk2
        midpoint_grid = copy(lattice.t_grid)
        for n in 1:(n_nodes - 1)
            midpoint_grid[n] = 0.5 * (
                lattice.t_grid[n] + lattice.t_grid[n + 1]
            )
        end
        manufactured_kernel_table(
            midpoint_grid,
            rates,
            time_origin = time_origin,
            time_slope = time_slope,
        )
    else
        nothing
    end

    solution = _evolve_stepwise_rg(
        lattice = lattice,
        boundary_values = boundary_values,
        kernel_table = kernel_table,
        midpoint_kernel_table = midpoint_kernel_table,
        method = method,
        store_history = true,
    )
    squared_error = (
        sum(abs2, solution.jq_history .- exact_q) +
        sum(abs2, solution.jg_history .- exact_g)
    )
    squared_reference = sum(abs2, exact_q) + sum(abs2, exact_g)
    maximum_error = max(
        maximum(abs, solution.jq_history .- exact_q),
        maximum(abs, solution.jg_history .- exact_g),
    )

    return (
        maximum_delta_t = maximum(diff(lattice.t_grid)),
        relative_l2 = sqrt(squared_error / squared_reference),
        maximum = maximum_error,
    )
end

@testset "Nonlocal matrix-exponential evolution" begin
    for method in (:euler, :rk2)
        coarse = manufactured_nonlocal_error(17, method)
        fine = manufactured_nonlocal_error(33, method)
        observed_order = log(coarse.relative_l2 / fine.relative_l2) /
                         log(coarse.maximum_delta_t / fine.maximum_delta_t)

        @test isfinite(coarse.relative_l2)
        @test isfinite(fine.relative_l2)
        @test fine.relative_l2 < coarse.relative_l2
        @test method == :euler ? observed_order > 0.9 : observed_order > 1.8
    end

    euler = manufactured_nonlocal_error(33, :euler)
    rk2 = manufactured_nonlocal_error(33, :rk2)
    @test rk2.relative_l2 < euler.relative_l2
    @test rk2.maximum < euler.maximum
end
