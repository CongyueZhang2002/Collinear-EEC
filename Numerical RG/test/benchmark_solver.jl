using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

const BENCHMARK_NF_SCHEME = :VFNS
n_nodes = isempty(ARGS) ? 1000 : parse(Int64, ARGS[1])
if length(ARGS) >= 2
    BLAS.set_num_threads(parse(Int64, ARGS[2]))
end

bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)

boundary_func(; b, bstar, mu_start) = (
    exp(-2.0 * b),
    0.5 * exp(-3.0 * b),
)

println("Julia version: ", VERSION)
println("Julia threads: ", Threads.nthreads())
println("BLAS threads: ", BLAS.get_num_threads())
println("Warmup started")
flush(stdout)

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
warmup_mu = sqrt(
    warmup_solution.lattice.mu_i_grid[1] *
    warmup_solution.lattice.mu_i_grid[end]
)
warmup_solution(0.5, warmup_mu)
warmup_solution = nothing
GC.gc()

println("Warmup finished")
println("Timed $n_nodes-node solve started")
flush(stdout)

timed_result = @timed solve_jet_rg(
    n_nodes = n_nodes,
    b_min = 0.001,
    b_max = 30.0,
    bstar_func = bstar_func,
    boundary_func = boundary_func,
    order = 2,
    nf_scheme = BENCHMARK_NF_SCHEME,
    method = :rk2,
)
solution = timed_result.value

probe_mu = sqrt(solution.lattice.mu_i_grid[1] * solution.lattice.mu_i_grid[end])
probe = solution(1.0, probe_mu)

@printf("RESULT elapsed_seconds=%.6f\n", timed_result.time)
@printf("RESULT gc_seconds=%.6f\n", timed_result.gctime)
@printf("RESULT allocated_gib=%.6f\n", timed_result.bytes / 1024.0^3)
@printf("RESULT history_shape=%dx%d\n", size(solution.jq_history)...)
@printf("RESULT mu_range=[%.12g, %.12g]\n",
    solution.lattice.mu_i_grid[1],
    solution.lattice.mu_i_grid[end],
)
@printf("RESULT probe_b=1 probe_mu=%.12g jq=%.12g jg=%.12g\n",
    probe_mu,
    probe.jq,
    probe.jg,
)
flush(stdout)
