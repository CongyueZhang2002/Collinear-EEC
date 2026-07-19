using QuadGK

function continuum_physical_rhs(
    b::Float64;
    a_s::Float64,
    order::Int64,
    nf::Int64,
    jq_func::Function,
    jg_func::Function,
)

    jq_current = jq_func(b)
    jg_current = jg_func(b)
    endpoint = timelike_splitting_convolution_func(
        y = 0.5, as = a_s, order = order, nf = nf,
    )

    function q_integrand(y)
        components = timelike_splitting_convolution_func(
            y = y, as = a_s, order = order, nf = nf,
        )
        jq_weighted = y^2 * jq_func(b / y)
        jg_weighted = y^2 * jg_func(b / y)
        return (
            jq_weighted * components.PqqReg +
            (
                jq_weighted * components.PqqD0 -
                jq_current * endpoint.PqqD0_at_1
            ) / (1.0 - y) +
            (
                jq_weighted * components.PqqD1 -
                jq_current * endpoint.PqqD1_at_1
            ) * log1p(-y) / (1.0 - y) +
            jg_weighted * components.Pgq
        )
    end

    function g_integrand(y)
        components = timelike_splitting_convolution_func(
            y = y, as = a_s, order = order, nf = nf,
        )
        jq_weighted = y^2 * jq_func(b / y)
        jg_weighted = y^2 * jg_func(b / y)
        return (
            jq_weighted * components.Pqg +
            jg_weighted * components.PggReg +
            (
                jg_weighted * components.PggD0 -
                jg_current * endpoint.PggD0_at_1
            ) / (1.0 - y) +
            (
                jg_weighted * components.PggD1 -
                jg_current * endpoint.PggD1_at_1
            ) * log1p(-y) / (1.0 - y)
        )
    end

    integration_points = (0.0, 0.5, 0.9, 0.99, 0.9999, 1.0)
    rhs_q = (
        jq_current * endpoint.PqqDelta_at_1 +
        quadgk(q_integrand, integration_points..., rtol = 2.0e-10)[1]
    )
    rhs_g = (
        jg_current * endpoint.PggDelta_at_1 +
        quadgk(g_integrand, integration_points..., rtol = 2.0e-10)[1]
    )
    return (q = rhs_q, g = rhs_g)
end

function physical_convolution_error(
    n_nodes::Int64,
    order::Int64,
    nf_scheme::Symbol,
)
    b_min = 1.0e-3
    b_max = 30.0
    ell_grid = collect(range(log(b_max), log(b_min), length = n_nodes))
    b_grid = exp.(ell_grid)
    delta_ell = (log(b_max) - log(b_min)) / (n_nodes - 1)
    mu = 20.0
    time_grid = collect(range(2.0 * log(mu), 2.0 * log(21.0), length = n_nodes))
    kernels = build_splitting_kernel_table(
        delta_ell = delta_ell,
        n_nodes = n_nodes,
        t_grid = time_grid,
        order = order,
        nf_scheme = nf_scheme,
    )

    jq_func(b) = exp(-2.0 * b) * (1.0 + 0.3 * b)
    jg_func(b) = 0.7 * exp(-1.5 * b) * (1.0 + 0.2 * b)
    jq = jq_func.(b_grid)
    jg = jg_func.(b_grid)
    node_index = (n_nodes + 1) ÷ 2
    lattice_rhs = _build_node_rhs(jq, jg, kernels, 1, node_index)
    a_s = alpha_s_func_numerical(
        mu_f = mu,
        order = 4,
        nf_scheme = nf_scheme,
    ) / (4.0 * pi)
    continuum_rhs = continuum_physical_rhs(
        b_grid[node_index],
        a_s = a_s,
        order = order,
        nf = nf_func(mu; scheme = nf_scheme),
        jq_func = jq_func,
        jg_func = jg_func,
    )

    q_relative_error = abs(lattice_rhs[1] - continuum_rhs.q) /
                       max(abs(continuum_rhs.q), eps())
    g_relative_error = abs(lattice_rhs[2] - continuum_rhs.g) /
                       max(abs(continuum_rhs.g), eps())
    return (
        delta_ell = delta_ell,
        q_relative_error = q_relative_error,
        g_relative_error = g_relative_error,
    )
end

@testset "Smooth physical convolution convergence" begin
    nf_scheme = :VFNS
    for order in 0:2
        coarse = physical_convolution_error(65, Int64(order), nf_scheme)
        medium = physical_convolution_error(129, Int64(order), nf_scheme)
        fine = physical_convolution_error(257, Int64(order), nf_scheme)

        @test medium.q_relative_error < coarse.q_relative_error
        @test fine.q_relative_error < medium.q_relative_error
        @test medium.g_relative_error < coarse.g_relative_error
        @test fine.g_relative_error < medium.g_relative_error

        q_order = log(medium.q_relative_error / fine.q_relative_error) /
                  log(medium.delta_ell / fine.delta_ell)
        g_order = log(medium.g_relative_error / fine.g_relative_error) /
                  log(medium.delta_ell / fine.delta_ell)
        @test q_order > 0.8
        @test g_order > 0.8
    end
end
