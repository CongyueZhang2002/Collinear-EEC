struct SplittingComponentGrid
    shift::Vector{Int64}
    y::Vector{Float64}
    weight::Vector{Float64}
    PqqReg::Vector{Float64}
    PqqDelta_at_1::Float64
    PqqD0::Vector{Float64}
    PqqD0_at_1::Float64
    PqqD1::Vector{Float64}
    PqqD1_at_1::Float64
    Pqg::Vector{Float64}
    Pgq::Vector{Float64}
    PggReg::Vector{Float64}
    PggDelta_at_1::Float64
    PggD0::Vector{Float64}
    PggD0_at_1::Float64
    PggD1::Vector{Float64}
    PggD1_at_1::Float64
    q_from_q_kernel::Vector{Float64}
    q_from_g_kernel::Vector{Float64}
    g_from_q_kernel::Vector{Float64}
    g_from_g_kernel::Vector{Float64}
    q_self_sub::Vector{Float64}
    g_self_sub::Vector{Float64}
end

function set_alpha_s_runtime!(; alpha_s_Z_value::Float64, alpha_s_loops::Int64)
    global αs_Z = alpha_s_Z_value
    global nloops_αs = alpha_s_loops
    return nothing
end

function build_splitting_component_grid(;
    delta_ell::Float64,
    n_nodes::Int64,
    mu::Float64,
    order::Int64,
    alpha_s_Z_value::Float64 = 0.118,
    alpha_s_loops::Int64 = 4,
)

    set_alpha_s_runtime!(alpha_s_Z_value = alpha_s_Z_value, alpha_s_loops = alpha_s_loops)

    a_s = alpha_s_func(mu) / (4 * pi)

    n_shift = max(n_nodes - 1, 0)
    shift = collect(1:n_shift)
    y = exp.(-delta_ell .* shift)
    weight = delta_ell .* y.^3

    # The *_at_1 coefficients are independent of y, so a safe y < 1 point is enough.
    y_ref = exp(-delta_ell)
    at_1_components = spliting_convolution_func(y = y_ref, as = a_s, order = order)

    PqqReg = zeros(Float64, n_shift)
    PqqD0 = zeros(Float64, n_shift)
    PqqD1 = zeros(Float64, n_shift)
    Pqg = zeros(Float64, n_shift)
    Pgq = zeros(Float64, n_shift)
    PggReg = zeros(Float64, n_shift)
    PggD0 = zeros(Float64, n_shift)
    PggD1 = zeros(Float64, n_shift)

    for k in 1:n_shift
        components = spliting_convolution_func(y = y[k], as = a_s, order = order)
        PqqReg[k] = components.PqqReg
        PqqD0[k] = components.PqqD0
        PqqD1[k] = components.PqqD1
        Pqg[k] = components.Pqg
        Pgq[k] = components.Pgq
        PggReg[k] = components.PggReg
        PggD0[k] = components.PggD0
        PggD1[k] = components.PggD1
    end

    one_minus_y = 1 .- y
    log_one_minus_y = log.(one_minus_y)

    q_from_q_kernel = weight .* (
        PqqReg .+
        PqqD0 ./ one_minus_y .+
        PqqD1 .* log_one_minus_y ./ one_minus_y
    )
    q_from_g_kernel = weight .* Pgq

    g_from_q_kernel = weight .* Pqg
    g_from_g_kernel = weight .* (
        PggReg .+
        PggD0 ./ one_minus_y .+
        PggD1 .* log_one_minus_y ./ one_minus_y
    )

    q_self_sub = cumsum(
        weight .* (
            at_1_components.PqqD0_at_1 ./ one_minus_y .+
            at_1_components.PqqD1_at_1 .* log_one_minus_y ./ one_minus_y
        )
    )
    g_self_sub = cumsum(
        weight .* (
            at_1_components.PggD0_at_1 ./ one_minus_y .+
            at_1_components.PggD1_at_1 .* log_one_minus_y ./ one_minus_y
        )
    )

    return SplittingComponentGrid(
        shift,
        y,
        weight,
        PqqReg,
        at_1_components.PqqDelta_at_1,
        PqqD0,
        at_1_components.PqqD0_at_1,
        PqqD1,
        at_1_components.PqqD1_at_1,
        Pqg,
        Pgq,
        PggReg,
        at_1_components.PggDelta_at_1,
        PggD0,
        at_1_components.PggD0_at_1,
        PggD1,
        at_1_components.PggD1_at_1,
        q_from_q_kernel,
        q_from_g_kernel,
        g_from_q_kernel,
        g_from_g_kernel,
        q_self_sub,
        g_self_sub,
    )
end
