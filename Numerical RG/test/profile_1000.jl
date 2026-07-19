using Profile

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)

boundary_func(; b, bstar, mu_start) = (
    exp(-2.0 * b),
    0.5 * exp(-3.0 * b),
)

warmup_solution = solve_jet_rg(
    n_nodes = 12,
    b_min = 0.1,
    b_max = 1.0,
    bstar_func = bstar_func,
    boundary_func = boundary_func,
    order = 2,
    method = :rk2,
)
warmup_solution(0.5, sqrt(
    warmup_solution.lattice.mu_i_grid[1] *
    warmup_solution.lattice.mu_i_grid[end]
))
warmup_solution = nothing
GC.gc()

Profile.clear()
Profile.init(n = 10^7, delay = 0.001)

solution = nothing
elapsed = @elapsed begin
    Profile.@profile solution = solve_jet_rg(
        n_nodes = 1000,
        b_min = 0.001,
        b_max = 30.0,
        bstar_func = bstar_func,
        boundary_func = boundary_func,
        order = 2,
        method = :rk2,
    )
end

println("PROFILE elapsed_seconds=", elapsed)
println("PROFILE samples=", length(Profile.fetch()))
Profile.print(format = :flat, sortedby = :count, mincount = 20)
