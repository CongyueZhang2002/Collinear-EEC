using Printf

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

function build_local_delta_table(
    n_nodes::Int64;
    q_rate::Float64,
    g_rate::Float64,
)

    n_shift = n_nodes - 1
    return SplittingKernelTable(
        fill(q_rate, n_nodes),
        fill(g_rate, n_nodes),
        zeros(Float64, n_nodes, n_shift),
        zeros(Float64, n_nodes, n_shift),
        zeros(Float64, n_nodes, n_shift),
        zeros(Float64, n_nodes, n_shift),
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
    )
end

function local_delta_error(
    n_nodes::Int64,
    method::Symbol;
    q_rate::Float64 = 0.4,
    g_rate::Float64 = -0.3,
)

    lattice = build_lattice_grid(
        n_nodes = n_nodes,
        b_min = 0.5,
        b_max = 2.0,
        bstar_func = identity,
    )
    boundary_func(; b, bstar, mu_start) = (1.0 + 0.2 * b, 0.7 + 0.1 * b)
    boundary_values = build_boundary_values(
        grid = lattice,
        boundary_func = boundary_func,
    )
    kernel_table = build_local_delta_table(
        n_nodes,
        q_rate = q_rate,
        g_rate = g_rate,
    )

    solution = _evolve_stepwise_rg(
        lattice = lattice,
        boundary_values = boundary_values,
        kernel_table = kernel_table,
        midpoint_kernel_table = method == :rk2 ? kernel_table : nothing,
        method = method,
        store_history = true,
    )

    squared_error = 0.0
    maximum_error = 0.0
    n_values = 0

    for n in 1:n_nodes
        for i in 1:n_nodes
            delta_t = lattice.t_grid[n] - lattice.t_grid[i]
            jq_exact = boundary_values.jq_bc[i] * exp(q_rate * delta_t)
            jg_exact = boundary_values.jg_bc[i] * exp(g_rate * delta_t)
            jq_error = abs(solution.jq_history[i, n] / jq_exact - 1.0)
            jg_error = abs(solution.jg_history[i, n] / jg_exact - 1.0)

            squared_error += jq_error^2 + jg_error^2
            maximum_error = max(maximum_error, jq_error, jg_error)
            n_values += 2
        end
    end

    return (
        delta_t = lattice.t_grid[2] - lattice.t_grid[1],
        rms = sqrt(squared_error / n_values),
        maximum = maximum_error,
    )
end

function triangle_relative_metrics(
    solution::StepwiseRGSolution,
    component::Symbol,
    region::Symbol,
    maximum_b::Float64,
)

    history = component == :q ? solution.jq_history : solution.jg_history
    squared_error = 0.0
    maximum_error = 0.0
    n_values = 0

    for n in 1:solution.lattice.n_nodes
        for i in 1:solution.lattice.n_nodes
            b = solution.lattice.b_grid[i]
            if b > maximum_b
                continue
            end
            if region == :upper && n < i
                continue
            end
            if region == :lower && n > i
                continue
            end

            relative_error = abs(history[i, n] / b - 1.0)
            squared_error += relative_error^2
            maximum_error = max(maximum_error, relative_error)
            n_values += 1
        end
    end

    return (
        rms = sqrt(squared_error / n_values),
        maximum = maximum_error,
        n_values = n_values,
    )
end

function run_momentum_test(
    n_nodes::Int64,
    order::Int64;
    grid_b_max::Float64 = 30.0,
)

    bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)
    momentum_boundary(; b, bstar, mu_start) = (b, b)

    solution = solve_jet_rg(
        n_nodes = n_nodes,
        b_min = 0.001,
        b_max = grid_b_max,
        bstar_func = bstar_func,
        boundary_func = momentum_boundary,
        order = order,
        method = :rk2,
    )

    for maximum_b in (0.03, 0.3)
        for region in (:upper, :lower)
            q_metrics = triangle_relative_metrics(solution, :q, region, maximum_b)
            g_metrics = triangle_relative_metrics(solution, :g, region, maximum_b)
            @printf(
                "MOMENTUM n_nodes=%d grid_b_max=%.3g order=%d b_max_eval=%.3g region=%s q_rms=%.8e q_max=%.8e g_rms=%.8e g_max=%.8e\n",
                n_nodes,
                grid_b_max,
                order,
                maximum_b,
                string(region),
                q_metrics.rms,
                q_metrics.maximum,
                g_metrics.rms,
                g_metrics.maximum,
            )
        end
    end
end

println("Julia threads: ", Threads.nthreads())

for method in (:euler, :rk2)
    previous = nothing
    for n_nodes in (17, 33, 65, 129)
        metrics = local_delta_error(n_nodes, method)
        observed_order = if previous === nothing
            NaN
        else
            log(previous.rms / metrics.rms) / log(previous.delta_t / metrics.delta_t)
        end
        @printf(
            "LOCAL method=%s n_nodes=%d dt=%.8e rms=%.8e max=%.8e observed_order=%.6f\n",
            string(method),
            n_nodes,
            metrics.delta_t,
            metrics.rms,
            metrics.maximum,
            observed_order,
        )
        previous = metrics
    end
end

momentum_cases = if isempty(ARGS)
    [(n_nodes = 1000, grid_b_max = 30.0)]
else
    map(ARGS) do argument
        fields = split(argument, ":")
        (
            n_nodes = parse(Int64, fields[1]),
            grid_b_max = length(fields) == 1 ? 30.0 : parse(Float64, fields[2]),
        )
    end
end

for case in momentum_cases
    run_momentum_test(
        case.n_nodes,
        0,
        grid_b_max = case.grid_b_max,
    )
    run_momentum_test(
        case.n_nodes,
        2,
        grid_b_max = case.grid_b_max,
    )
end
