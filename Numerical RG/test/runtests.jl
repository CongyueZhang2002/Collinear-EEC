using Test

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))
include(joinpath(@__DIR__, "timelike_splitting_tests.jl"))
include(joinpath(@__DIR__, "strong_coupling_tests.jl"))

function test_components(actual, expected)
    @test propertynames(actual) == propertynames(expected)
    for name in propertynames(actual)
        @test getproperty(actual, name) ≈ getproperty(expected, name)
    end
end

@testset "Lattice grid" begin
    bstar_func(b) = b / sqrt(1.0 + b^2)
    grid = build_lattice_grid(
        n_nodes = 7,
        b_min = 0.1,
        b_max = 1.0,
        bstar_func = bstar_func,
    )

    @test grid.n_nodes == 7
    @test axes(grid.ell_grid, 1) == Base.OneTo(grid.n_nodes)
    @test grid.b_grid ≈ exp.(grid.ell_grid)
    @test grid.bstar_grid ≈ bstar_func.(grid.b_grid)
    @test grid.mu_i_grid ≈ b0 ./ grid.bstar_grid
    @test grid.t_grid ≈ log.(grid.mu_i_grid.^2)
    @test all(<(0.0), diff(grid.ell_grid))
    @test all(<(0.0), diff(grid.b_grid))
    @test all(>(0.0), diff(grid.t_grid))

    notebook_grid = build_lattice_grid(
        n_nodes = 1000,
        b_min = 5.0e-4,
        b_max = 30.0,
        bstar_func = b -> b / sqrt(1.0 + (b / 1.122918)^2),
    )
    @test axes(notebook_grid.t_grid, 1) == Base.OneTo(notebook_grid.n_nodes)
    @test all(>(0.0), diff(notebook_grid.t_grid))

    @test_throws ArgumentError build_lattice_grid(
        n_nodes = 1,
        b_min = 0.1,
        b_max = 1.0,
        bstar_func = identity,
    )
    @test_throws ArgumentError build_lattice_grid(
        n_nodes = 3,
        b_min = 1.0,
        b_max = 1.0,
        bstar_func = identity,
    )
    @test_throws ArgumentError build_lattice_grid(
        n_nodes = 3,
        b_min = 0.1,
        b_max = 1.0,
        bstar_func = b -> 1.0 / b,
    )
    @test_throws DomainError build_lattice_grid(
        n_nodes = 3,
        b_min = 0.1,
        b_max = 1.0,
        bstar_func = b -> 0.0,
    )
    @test_throws DomainError build_lattice_grid(
        n_nodes = 3,
        b_min = 0.1,
        b_max = 1.0,
        bstar_func = b -> NaN,
    )

    @test_throws DomainError build_boundary_values(
        grid = grid,
        boundary_func = (; b, bstar, mu_start) -> (NaN, 1.0),
    )
end

@testset "Perturbative splitting series" begin
    y = 0.37
    a_s = 0.013
    nf = nf_func(10.0; scheme = :VFNS)

    Pqq0Reg, Pqq0Delta, Pqq0D0 = timelike_splittings_func(
        x = y, type = "Pqq0", nf = nf,
    )
    Pqg0 = timelike_splittings_func(x = y, type = "Pqg0", nf = nf)
    Pgq0 = timelike_splittings_func(x = y, type = "Pgq0", nf = nf)
    Pgg0Reg, Pgg0Delta, Pgg0D0 = timelike_splittings_func(
        x = y, type = "Pgg0", nf = nf,
    )

    expected_lo = (
        PqqReg = a_s * Pqq0Reg,
        PqqDelta_at_1 = a_s * Pqq0Delta,
        PqqD0 = a_s * Pqq0D0,
        PqqD0_at_1 = a_s * Pqq0D0,
        PqqD1 = 0.0,
        PqqD1_at_1 = 0.0,
        Pqg = a_s * Pqg0,
        Pgq = a_s * Pgq0,
        PggReg = a_s * Pgg0Reg,
        PggDelta_at_1 = a_s * Pgg0Delta,
        PggD0 = a_s * Pgg0D0,
        PggD0_at_1 = a_s * Pgg0D0,
        PggD1 = 0.0,
        PggD1_at_1 = 0.0,
    )
    lo = timelike_splitting_convolution_func(y = y, as = a_s, order = 0, nf = nf)
    test_components(lo, expected_lo)

    Pqq1Reg, Pqq1Delta, Pqq1D0, Pqq1D1 = timelike_splittings_func(
        x = y, type = "Pqq1", nf = nf,
    )
    Pqg1 = timelike_splittings_func(x = y, type = "Pqg1", nf = nf)
    Pgq1 = timelike_splittings_func(x = y, type = "Pgq1", nf = nf)
    Pgg1Reg, Pgg1Delta, Pgg1D0, Pgg1D1 = timelike_splittings_func(
        x = y, type = "Pgg1", nf = nf,
    )
    Pqq1D0At1 = (-16.0 / 27.0) * (-201.0 + 10.0 * nf + 9.0 * pi^2)
    Pgg1D0At1 = (-4.0 / 3.0) * (-201.0 + 10.0 * nf + 9.0 * pi^2)
    a_s2 = a_s^2

    expected_nlo = (
        PqqReg = lo.PqqReg + a_s2 * Pqq1Reg,
        PqqDelta_at_1 = lo.PqqDelta_at_1 + a_s2 * Pqq1Delta,
        PqqD0 = lo.PqqD0 + a_s2 * Pqq1D0,
        PqqD0_at_1 = lo.PqqD0_at_1 + a_s2 * Pqq1D0At1,
        PqqD1 = a_s2 * Pqq1D1,
        PqqD1_at_1 = 0.0,
        Pqg = lo.Pqg + a_s2 * Pqg1,
        Pgq = lo.Pgq + a_s2 * Pgq1,
        PggReg = lo.PggReg + a_s2 * Pgg1Reg,
        PggDelta_at_1 = lo.PggDelta_at_1 + a_s2 * Pgg1Delta,
        PggD0 = lo.PggD0 + a_s2 * Pgg1D0,
        PggD0_at_1 = lo.PggD0_at_1 + a_s2 * Pgg1D0At1,
        PggD1 = a_s2 * Pgg1D1,
        PggD1_at_1 = 0.0,
    )
    nlo = timelike_splitting_convolution_func(y = y, as = a_s, order = 1, nf = nf)
    test_components(nlo, expected_nlo)

    Pqq2Reg, Pqq2Delta, Pqq2D0 = timelike_splittings_func(
        x = y, type = "Pqq2", nf = nf,
    )
    Pqg2 = timelike_splittings_func(x = y, type = "Pqg2", nf = nf)
    Pgq2 = timelike_splittings_func(x = y, type = "Pgq2", nf = nf)
    Pgg2Reg, Pgg2Delta, Pgg2D0 = timelike_splittings_func(
        x = y, type = "Pgg2", nf = nf,
    )
    a_s3 = a_s^3

    expected_nnlo = (
        PqqReg = nlo.PqqReg + a_s3 * Pqq2Reg,
        PqqDelta_at_1 = nlo.PqqDelta_at_1 + a_s3 * Pqq2Delta,
        PqqD0 = nlo.PqqD0 + a_s3 * Pqq2D0,
        PqqD0_at_1 = nlo.PqqD0_at_1 + a_s3 * Pqq2D0,
        PqqD1 = nlo.PqqD1,
        PqqD1_at_1 = nlo.PqqD1_at_1,
        Pqg = nlo.Pqg + a_s3 * Pqg2,
        Pgq = nlo.Pgq + a_s3 * Pgq2,
        PggReg = nlo.PggReg + a_s3 * Pgg2Reg,
        PggDelta_at_1 = nlo.PggDelta_at_1 + a_s3 * Pgg2Delta,
        PggD0 = nlo.PggD0 + a_s3 * Pgg2D0,
        PggD0_at_1 = nlo.PggD0_at_1 + a_s3 * Pgg2D0,
        PggD1 = nlo.PggD1,
        PggD1_at_1 = nlo.PggD1_at_1,
    )
    nnlo = timelike_splitting_convolution_func(y = y, as = a_s, order = 2, nf = nf)
    test_components(nnlo, expected_nnlo)

    @test_throws DomainError timelike_splitting_convolution_func(
        y = 0.0, as = a_s, order = 0, nf = nf,
    )
    @test_throws DomainError timelike_splitting_convolution_func(
        y = 1.0, as = a_s, order = 0, nf = nf,
    )
    @test_throws DomainError timelike_splitting_convolution_func(
        y = y, as = -a_s, order = 0, nf = nf,
    )
    @test_throws DomainError timelike_splitting_convolution_func(
        y = y, as = NaN, order = 0, nf = nf,
    )
    @test_throws ArgumentError timelike_splitting_convolution_func(
        y = y, as = a_s, order = 3, nf = nf,
    )
end

include("paper_moment_tests.jl")

@testset "Precomputed splitting table" begin
    nf = nf_func(10.0; scheme = :nf5)
    grid = build_lattice_grid(
        n_nodes = 5,
        b_min = 0.2,
        b_max = 1.0,
        bstar_func = identity,
    )

    for order in 0:2
        table = build_splitting_kernel_table(
            delta_ell = grid.delta_ell,
            n_nodes = grid.n_nodes,
            t_grid = grid.t_grid,
            order = order,
            nf_scheme = :nf5,
        )

        @test size(table.qq) == (grid.n_nodes, grid.n_nodes - 1)
        @test size(table.qg) == size(table.qq)
        @test size(table.gq) == size(table.qq)
        @test size(table.gg) == size(table.qq)
        @test axes(table.qq_d0_at_1_vec, 1) == Base.OneTo(grid.n_nodes)
        @test axes(table.qq_d1_at_1_vec, 1) == Base.OneTo(grid.n_nodes)
        @test axes(table.gg_d0_at_1_vec, 1) == Base.OneTo(grid.n_nodes)
        @test axes(table.gg_d1_at_1_vec, 1) == Base.OneTo(grid.n_nodes)
        @test axes(table.d0_sub_weight_vec, 1) == Base.OneTo(grid.n_nodes - 1)
        @test axes(table.d1_sub_weight_vec, 1) == Base.OneTo(grid.n_nodes - 1)
        @test axes(table.d0_closure_factor_vec, 1) == Base.OneTo(grid.n_nodes)
        @test axes(table.d1_closure_factor_vec, 1) == Base.OneTo(grid.n_nodes)
        @test all(isfinite, table.qq)
        @test all(isfinite, table.qg)
        @test all(isfinite, table.gq)
        @test all(isfinite, table.gg)
        @test all(isfinite, table.qq_d0_at_1_vec)
        @test all(isfinite, table.qq_d1_at_1_vec)
        @test all(isfinite, table.gg_d0_at_1_vec)
        @test all(isfinite, table.gg_d1_at_1_vec)
    end

    scheme_scales = [1.1, 2.0, 5.0, 100.0, 300.0]
    scheme_t_grid = log.(scheme_scales.^2)
    for nf_scheme in NF_SCHEMES
        scheme_table = build_splitting_kernel_table(
            delta_ell = grid.delta_ell,
            n_nodes = grid.n_nodes,
            t_grid = scheme_t_grid,
            order = 2,
            nf_scheme = nf_scheme,
        )
        @test all(isfinite, scheme_table.qq)
        @test all(isfinite, scheme_table.qg)
        @test all(isfinite, scheme_table.gq)
        @test all(isfinite, scheme_table.gg)
    end

    vfns_table = build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = scheme_t_grid,
        order = 0,
        nf_scheme = :VFNS,
    )
    for n in eachindex(scheme_scales)
        mu = scheme_scales[n]
        alpha_s = alpha_s_func_numerical(
            mu_f = mu,
            order = 4,
            nf_scheme = :VFNS,
        )
        expected = timelike_splitting_convolution_func(
            y = exp(-grid.delta_ell),
            as = alpha_s / (4.0 * pi),
            order = 0,
            nf = nf_func(mu; scheme = :VFNS),
        )
        @test vfns_table.qq_delta_at_1_vec[n] ≈ expected.PqqDelta_at_1
        @test vfns_table.gg_delta_at_1_vec[n] ≈ expected.PggDelta_at_1
    end

    table = build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 1,
        nf_scheme = :nf5,
    )
    y_vec = exp.(-grid.delta_ell .* collect(1:(grid.n_nodes - 1)))
    one_minus_y_vec = 1.0 .- y_vec
    tail_log = log1p(-y_vec[end])
    d0_total_factor = sum(grid.delta_ell .* y_vec ./ one_minus_y_vec) - tail_log
    d1_total_factor = (
        sum(grid.delta_ell .* y_vec .* log.(one_minus_y_vec) ./ one_minus_y_vec) -
        0.5 * tail_log^2
    )
    alpha_s = alpha_s_func_numerical(
        mu_f = grid.mu_i_grid[1],
        order = 4,
        nf_scheme = :nf5,
    )
    at_1 = timelike_splitting_convolution_func(
        y = exp(-grid.delta_ell),
        as = alpha_s / (4.0 * pi),
        order = 1,
        nf = nf,
    )
    @test table.qq_d0_at_1_vec[1] ≈ at_1.PqqD0_at_1
    @test table.qq_d1_at_1_vec[1] ≈ at_1.PqqD1_at_1
    @test table.gg_d0_at_1_vec[1] ≈ at_1.PggD0_at_1
    @test table.gg_d1_at_1_vec[1] ≈ at_1.PggD1_at_1
    @test table.d0_closure_factor_vec[1] ≈ d0_total_factor
    @test table.d1_closure_factor_vec[1] ≈ d1_total_factor
    @test table.d0_closure_factor_vec[end] ≈ -tail_log
    @test table.d1_closure_factor_vec[end] ≈ -0.5 * tail_log^2

    first_components = timelike_splitting_convolution_func(
        y = y_vec[1],
        as = alpha_s / (4.0 * pi),
        order = 1,
        nf = nf,
    )
    shifted_weight = grid.delta_ell * y_vec[1]^3
    @test table.qg[1, 1] ≈ shifted_weight * first_components.Pqg
    @test table.gq[1, 1] ≈ shifted_weight * first_components.Pgq

    @test_throws DimensionMismatch build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid[1:end-1],
        order = 0,
        nf_scheme = :nf5,
    )
    @test_throws ArgumentError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = reverse(grid.t_grid),
        order = 0,
        nf_scheme = :nf5,
    )
    @test_throws ArgumentError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 0,
        nf_scheme = :invalid,
    )
    @test_throws ArgumentError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 0,
        nf_scheme = :nf5,
        alpha_s_order = 0,
    )
    @test_throws ArgumentError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 0,
        nf_scheme = :nf5,
        alpha_s_order = 6,
    )
    @test_throws ArgumentError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 0,
        nf_scheme = :nf5,
        alpha_s_ref = -0.118,
    )
    @test_throws DomainError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = grid.n_nodes,
        t_grid = grid.t_grid,
        order = 0,
        nf_scheme = :nf5,
        alpha_s_max = 0.01,
    )
    @test_throws DomainError build_splitting_kernel_table(
        delta_ell = grid.delta_ell,
        n_nodes = 3,
        t_grid = log.([0.05, 0.06, 0.07].^2),
        order = 0,
        nf_scheme = :nf5,
        alpha_s_max = 1.0e8,
    )
end

@testset "Convolution channel and index mapping" begin
    kernels = SplittingKernelTable(
        [10.0],
        [20.0],
        reshape([1.0, 2.0], 1, 2),
        reshape([3.0, 4.0], 1, 2),
        reshape([5.0, 6.0], 1, 2),
        reshape([7.0, 8.0], 1, 2),
        [0.5],
        [0.0],
        [0.25],
        [0.0],
        [2.0, 3.0],
        [0.0, 0.0],
        [7.0, 11.0, 13.0],
        [0.0, 0.0, 0.0],
    )
    jq = [1.0, 2.0, 3.0]
    jg = [4.0, 5.0, 6.0]

    rhs = build_rhs(
        jq = jq,
        jg = jg,
        active = trues(3),
        kernels = kernels,
        time_index = 1,
    )
    @test rhs.rhs_q ≈ [6.5, 28.0, 56.0]
    @test rhs.rhs_g ≈ [73.0, 114.75, 170.0]

    inactive_rhs = build_rhs(
        jq = jq,
        jg = jg,
        active = falses(3),
        kernels = kernels,
        time_index = 1,
    )
    @test inactive_rhs.rhs_q == zeros(3)
    @test inactive_rhs.rhs_g == zeros(3)

    @test_throws DimensionMismatch build_rhs(
        jq = jq[1:2],
        jg = jg,
        active = trues(3),
        kernels = kernels,
        time_index = 1,
    )
    @test_throws BoundsError build_rhs(
        jq = jq,
        jg = jg,
        active = trues(3),
        kernels = kernels,
        time_index = 2,
    )
end

include("manufactured_nonlocal_tests.jl")
include("physical_convolution_tests.jl")

@testset "Stepwise solver" begin
    grid = build_lattice_grid(
        n_nodes = 5,
        b_min = 0.25,
        b_max = 1.0,
        bstar_func = identity,
    )
    boundary_func(; b, bstar, mu_start) = (b + bstar, mu_start)

    for method in (:euler, :rk2)
        result = solve_stepwise_rg(
            lattice = grid,
            boundary_func = boundary_func,
            order = 0,
            nf_scheme = :nf5,
            method = method,
            closure_check = :ignore,
        )

        @test all(result.active)
        @test all(isfinite, result.jq)
        @test all(isfinite, result.jg)
        @test size(result.jq_history) == (grid.n_nodes, grid.n_nodes)
        @test size(result.jg_history) == (grid.n_nodes, grid.n_nodes)
        @test all(isfinite, result.jq_history)
        @test all(isfinite, result.jg_history)
        for i in 1:grid.n_nodes
            @test result.jq_history[i, i] == result.jq_bc[i]
            @test result.jg_history[i, i] == result.jg_bc[i]
        end

        if method == :euler
            kernels = build_splitting_kernel_table(
                delta_ell = grid.delta_ell,
                n_nodes = grid.n_nodes,
                t_grid = result.time_grid,
                order = 0,
                nf_scheme = :nf5,
            )
            jq_state = zeros(Float64, grid.n_nodes)
            jg_state = zeros(Float64, grid.n_nodes)

            for i in 2:grid.n_nodes
                for n in 1:(i - 1)
                    jq_state[1:i] .= result.jq_history[1:i, n + 1]
                    jg_state[1:i] .= result.jg_history[1:i, n + 1]
                    rhs_q, rhs_g = _build_node_rhs(
                        jq_state,
                        jg_state,
                        kernels,
                        n + 1,
                        i,
                    )
                    dt = result.time_grid[n + 1] - result.time_grid[n]

                    @test result.jq_history[i, n] ≈ (
                        result.jq_history[i, n + 1] -
                        dt * rhs_q
                    )
                    @test result.jg_history[i, n] ≈ (
                        result.jg_history[i, n + 1] -
                        dt * rhs_g
                    )
                end
            end
        end
    end

    no_history = solve_stepwise_rg(
        lattice = grid,
        boundary_func = boundary_func,
        order = 0,
        nf_scheme = :nf5,
        method = :euler,
        store_history = false,
        closure_check = :ignore,
    )
    @test size(no_history.jq_history) == (0, 0)
    @test size(no_history.jg_history) == (0, 0)
    @test_throws ArgumentError no_history(
        grid.b_grid[1],
        grid.mu_i_grid[end],
    )

    @test_throws ArgumentError solve_stepwise_rg(
        lattice = grid,
        boundary_func = boundary_func,
        order = 0,
        nf_scheme = :nf5,
        method = :unknown,
    )
    @test_throws DomainError solve_stepwise_rg(
        lattice = grid,
        boundary_func = boundary_func,
        order = 0,
        nf_scheme = :nf5,
        alpha_s_max = 0.01,
        closure_check = :ignore,
    )
    @test_throws UndefKeywordError solve_stepwise_rg(
        lattice = grid,
        order = 0,
        nf_scheme = :nf5,
    )
    @test_throws UndefKeywordError solve_jet_rg(
        n_nodes = 3,
        b_min = 0.5,
        b_max = 1.0,
        bstar_func = identity,
        order = 0,
        nf_scheme = :nf5,
    )
    default_scheme_result = solve_jet_rg(
        n_nodes = 3,
        b_min = 0.5,
        b_max = 1.0,
        bstar_func = identity,
        boundary_func = boundary_func,
        order = 0,
        closure_check = :ignore,
    )
    @test all(isfinite, default_scheme_result.jq_history)
    @test_throws ArgumentError solve_jet_rg(
        n_nodes = 3,
        b_min = 0.5,
        b_max = 1.0,
        bstar_func = identity,
        boundary_func = boundary_func,
        order = 0,
        nf_scheme = :invalid,
        closure_check = :ignore,
    )
end

@testset "Zero-closure boundary checks" begin
    zero_boundary = (jq_bc = zeros(3), jg_bc = zeros(3))
    unsafe_boundary = (jq_bc = [0.1, 0.5, 1.0], jg_bc = [0.2, 0.4, 0.8])

    @test_logs check_zero_closure(
        zero_boundary,
        closure_check = :warn,
        closure_rtol = 1.0e-6,
        closure_atol = 1.0e-12,
    )
    @test_logs (:warn, r"Zero closure requires") check_zero_closure(
        unsafe_boundary,
        closure_check = :warn,
        closure_rtol = 1.0e-6,
        closure_atol = 1.0e-12,
    )
    @test_throws ArgumentError check_zero_closure(
        unsafe_boundary,
        closure_check = :error,
        closure_rtol = 1.0e-6,
        closure_atol = 1.0e-12,
    )
    @test check_zero_closure(
        unsafe_boundary,
        closure_check = :ignore,
        closure_rtol = 0.0,
        closure_atol = 0.0,
    ) === nothing
    @test_throws ArgumentError check_zero_closure(
        unsafe_boundary,
        closure_check = :invalid,
        closure_rtol = 1.0e-6,
        closure_atol = 1.0e-12,
    )
    @test_throws ArgumentError check_zero_closure(
        unsafe_boundary,
        closure_check = :warn,
        closure_rtol = -1.0,
        closure_atol = 0.0,
    )
end

@testset "Callable RG interpolation" begin
    n_nodes = 5
    b_min = 0.25
    b_max = 1.0
    boundary_func(; b, bstar, mu_start) = (1.0 + b, 2.0 + bstar)

    solution = solve_jet_rg(
        n_nodes = n_nodes,
        b_min = b_min,
        b_max = b_max,
        bstar_func = identity,
        boundary_func = boundary_func,
        order = 0,
        nf_scheme = :nf5,
        method = :rk2,
        closure_check = :ignore,
    )

    mu_min = solution.lattice.mu_i_grid[1]
    mu_max = solution.lattice.mu_i_grid[end]
    @test solution isa StepwiseRGSolution
    @test solution.time_grid == solution.lattice.t_grid
    @test solution.time_grid[1] ≈ 2.0 * log(mu_min)
    @test solution.time_grid[end] ≈ 2.0 * log(mu_max)
    @test size(solution.jq_history) == (n_nodes, n_nodes)
    @test size(solution.jg_history) == size(solution.jq_history)
    @test solution.jq == solution.jq_history[:, end]
    @test solution.jg == solution.jg_history[:, end]
    @test all(isfinite, solution.jq_history)
    @test all(isfinite, solution.jg_history)

    for i in 1:n_nodes
        @test mu_start(solution, solution.lattice.b_grid[i]) ≈ solution.lattice.mu_i_grid[i]
        for n in 1:size(solution.time_grid, 1)
            mu = exp(0.5 * solution.time_grid[n])
            value = solution(solution.lattice.b_grid[i], mu)
            @test value.jq ≈ solution.jq_history[i, n]
            @test value.jg ≈ solution.jg_history[i, n]
        end
    end

    b_index = 2
    b_fraction = 0.35
    ell = (
        (1.0 - b_fraction) * solution.lattice.ell_grid[b_index] +
        b_fraction * solution.lattice.ell_grid[b_index + 1]
    )
    b = exp(ell)

    t_boundary = (
        (1.0 - b_fraction) * solution.lattice.t_grid[b_index] +
        b_fraction * solution.lattice.t_grid[b_index + 1]
    )
    boundary_value = solution(b, exp(0.5 * t_boundary))
    @test boundary_value.jq ≈ (
        (1.0 - b_fraction) * solution.jq_history[b_index, b_index] +
        b_fraction * solution.jq_history[b_index + 1, b_index + 1]
    )
    @test boundary_value.jg ≈ (
        (1.0 - b_fraction) * solution.jg_history[b_index, b_index] +
        b_fraction * solution.jg_history[b_index + 1, b_index + 1]
    )

    lower_time_fraction = 0.2
    lower_t = (
        (1.0 - lower_time_fraction) * solution.time_grid[b_index] +
        lower_time_fraction * solution.time_grid[b_index + 1]
    )
    lower_value = solution(b, exp(0.5 * lower_t))
    @test lower_value.jq ≈ (
        (1.0 - b_fraction) * solution.jq_history[b_index, b_index] +
        (b_fraction - lower_time_fraction) *
        solution.jq_history[b_index + 1, b_index] +
        lower_time_fraction * solution.jq_history[b_index + 1, b_index + 1]
    )
    @test lower_value.jg ≈ (
        (1.0 - b_fraction) * solution.jg_history[b_index, b_index] +
        (b_fraction - lower_time_fraction) *
        solution.jg_history[b_index + 1, b_index] +
        lower_time_fraction * solution.jg_history[b_index + 1, b_index + 1]
    )

    diagonal_time_fraction = 0.7
    diagonal_t = (
        (1.0 - diagonal_time_fraction) * solution.time_grid[b_index] +
        diagonal_time_fraction * solution.time_grid[b_index + 1]
    )
    diagonal_value = solution(b, exp(0.5 * diagonal_t))
    @test diagonal_value.jq ≈ (
        (1.0 - diagonal_time_fraction) * solution.jq_history[b_index, b_index] +
        (diagonal_time_fraction - b_fraction) *
        solution.jq_history[b_index, b_index + 1] +
        b_fraction * solution.jq_history[b_index + 1, b_index + 1]
    )
    @test diagonal_value.jg ≈ (
        (1.0 - diagonal_time_fraction) * solution.jg_history[b_index, b_index] +
        (diagonal_time_fraction - b_fraction) *
        solution.jg_history[b_index, b_index + 1] +
        b_fraction * solution.jg_history[b_index + 1, b_index + 1]
    )

    time_index = b_index + 1
    time_fraction = 0.4
    interior_t = (
        (1.0 - time_fraction) * solution.time_grid[time_index] +
        time_fraction * solution.time_grid[time_index + 1]
    )
    interior_value = solution(b, exp(0.5 * interior_t))
    expected_q = (
        (1.0 - b_fraction) * (1.0 - time_fraction) *
        solution.jq_history[b_index, time_index] +
        b_fraction * (1.0 - time_fraction) *
        solution.jq_history[b_index + 1, time_index] +
        (1.0 - b_fraction) * time_fraction *
        solution.jq_history[b_index, time_index + 1] +
        b_fraction * time_fraction *
        solution.jq_history[b_index + 1, time_index + 1]
    )
    expected_g = (
        (1.0 - b_fraction) * (1.0 - time_fraction) *
        solution.jg_history[b_index, time_index] +
        b_fraction * (1.0 - time_fraction) *
        solution.jg_history[b_index + 1, time_index] +
        (1.0 - b_fraction) * time_fraction *
        solution.jg_history[b_index, time_index + 1] +
        b_fraction * time_fraction *
        solution.jg_history[b_index + 1, time_index + 1]
    )
    @test interior_value.jq ≈ expected_q
    @test interior_value.jg ≈ expected_g
    @test solution(b, exp(0.5 * interior_t); component = :q) ≈ expected_q
    @test solution(b, exp(0.5 * interior_t); component = :g) ≈ expected_g

    lower_time_index = 1
    lower_rectangle_fraction = 0.4
    lower_rectangle_t = (
        (1.0 - lower_rectangle_fraction) * solution.time_grid[lower_time_index] +
        lower_rectangle_fraction * solution.time_grid[lower_time_index + 1]
    )
    lower_rectangle_value = solution(b, exp(0.5 * lower_rectangle_t))
    lower_expected_q = (
        (1.0 - b_fraction) * (1.0 - lower_rectangle_fraction) *
        solution.jq_history[b_index, lower_time_index] +
        b_fraction * (1.0 - lower_rectangle_fraction) *
        solution.jq_history[b_index + 1, lower_time_index] +
        (1.0 - b_fraction) * lower_rectangle_fraction *
        solution.jq_history[b_index, lower_time_index + 1] +
        b_fraction * lower_rectangle_fraction *
        solution.jq_history[b_index + 1, lower_time_index + 1]
    )
    lower_expected_g = (
        (1.0 - b_fraction) * (1.0 - lower_rectangle_fraction) *
        solution.jg_history[b_index, lower_time_index] +
        b_fraction * (1.0 - lower_rectangle_fraction) *
        solution.jg_history[b_index + 1, lower_time_index] +
        (1.0 - b_fraction) * lower_rectangle_fraction *
        solution.jg_history[b_index, lower_time_index + 1] +
        b_fraction * lower_rectangle_fraction *
        solution.jg_history[b_index + 1, lower_time_index + 1]
    )
    @test lower_rectangle_value.jq ≈ lower_expected_q
    @test lower_rectangle_value.jg ≈ lower_expected_g

    below_boundary_t = 0.5 * (solution.time_grid[1] + t_boundary)
    below_boundary_value = solution(b, exp(0.5 * below_boundary_t))
    @test isfinite(below_boundary_value.jq)
    @test isfinite(below_boundary_value.jg)

    @test_throws DomainError solution(0.9 * b_min, mu_max)
    @test_throws DomainError solution(
        b,
        0.99 * exp(0.5 * solution.time_grid[1]),
    )
    @test_throws DomainError solution(b, 1.01 * mu_max)
    @test_throws ArgumentError solution(b, mu_max; component = :invalid)
end
