using Printf

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

const BENCHMARK_NF_SCHEME = :VFNS
n_nodes = isempty(ARGS) ? 1000 : parse(Int64, ARGS[1])

bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)

boundary_func(; b, bstar, mu_start) = (
    exp(-2.0 * b),
    0.5 * exp(-3.0 * b),
)

function build_rk2_kernel_tables(lattice::LatticeGrid)
    kernel_table = build_splitting_kernel_table(
        delta_ell = lattice.delta_ell,
        n_nodes = lattice.n_nodes,
        t_grid = lattice.t_grid,
        order = 2,
        nf_scheme = BENCHMARK_NF_SCHEME,
    )

    midpoint_grid = copy(lattice.t_grid)
    for n in 1:(lattice.n_nodes - 1)
        midpoint_grid[n] = 0.5 * (lattice.t_grid[n] + lattice.t_grid[n + 1])
    end
    midpoint_grid[end] = lattice.t_grid[end]

    midpoint_kernel_table = build_splitting_kernel_table(
        delta_ell = lattice.delta_ell,
        n_nodes = lattice.n_nodes,
        t_grid = midpoint_grid,
        order = 2,
        nf_scheme = BENCHMARK_NF_SCHEME,
    )

    return kernel_table, midpoint_kernel_table
end

# Compile every measured path on a small lattice first.
warmup_solution = solve_jet_rg(
    n_nodes = 12,
    b_min = 0.1,
    b_max = 1.0,
    bstar_func = bstar_func,
    boundary_func = boundary_func,
    order = 2,
    nf_scheme = BENCHMARK_NF_SCHEME,
    method = :rk2,
    closure_check = :ignore,
)
warmup_tables = build_rk2_kernel_tables(warmup_solution.lattice)
fill_lower_triangle!(
    copy(warmup_solution.jq_history),
    copy(warmup_solution.jg_history),
    warmup_solution.time_grid,
    warmup_tables...,
    :rk2,
    warmup_solution.lattice.n_nodes,
)
warmup_solution = nothing
warmup_tables = nothing
GC.gc()

total_result = @timed solve_jet_rg(
    n_nodes = n_nodes,
    b_min = 0.001,
    b_max = 30.0,
    bstar_func = bstar_func,
    boundary_func = boundary_func,
    order = 2,
    nf_scheme = BENCHMARK_NF_SCHEME,
    method = :rk2,
)
solution = total_result.value

lattice = solution.lattice
GC.gc()
kernel_result = @timed build_rk2_kernel_tables(lattice)
kernel_table, midpoint_kernel_table = kernel_result.value

jq_history = copy(solution.jq_history)
jg_history = copy(solution.jg_history)
GC.gc()
lower_result = @timed fill_lower_triangle!(
    jq_history,
    jg_history,
    solution.time_grid,
    kernel_table,
    midpoint_kernel_table,
    :rk2,
    n_nodes,
)

upper_and_setup = total_result.time - kernel_result.time - lower_result.time

println("Julia threads: ", Threads.nthreads())
@printf("COMPONENT total_seconds=%.6f\n", total_result.time)
@printf("COMPONENT kernels_seconds=%.6f\n", kernel_result.time)
@printf("COMPONENT lower_seconds=%.6f\n", lower_result.time)
@printf("COMPONENT upper_and_setup_seconds=%.6f\n", upper_and_setup)
