function unit_boundary_value(; b::Float64, bstar::Float64, mu_start::Float64, t_start::Float64)
    return 1.0, 1.0
end

function activate_nodes!(
    jq::Vector{Float64},
    jg::Vector{Float64},
    active::AbstractVector{Bool},
    jq_bc::AbstractVector{Float64},
    jg_bc::AbstractVector{Float64},
    t_start::AbstractVector{Float64},
    t::Float64;
    activation_tol::Float64 = 1e-12,
)

    for i in eachindex(active)
        if !active[i] && t + activation_tol >= t_start[i]
            jq[i] = jq_bc[i]
            jg[i] = jg_bc[i]
            active[i] = true
        end
    end

    return nothing
end

function solve_stepwise_rg(;
    lattice::LatticeGrid,
    time_grid::AbstractVector{Float64},
    boundary_func::Function = unit_boundary_value,
    order::Int64,
    alpha_s_Z_value::Float64 = 0.118,
    alpha_s_loops::Int64 = 4,
    method::Symbol = :euler,
    store_history::Bool = true,
)

    boundary_values = build_boundary_values(grid = lattice, boundary_func = boundary_func)

    jq = zeros(Float64, lattice.n_nodes)
    jg = zeros(Float64, lattice.n_nodes)
    active = falses(lattice.n_nodes)

    if store_history
        jq_history = zeros(Float64, lattice.n_nodes, length(time_grid))
        jg_history = zeros(Float64, lattice.n_nodes, length(time_grid))
    else
        jq_history = zeros(Float64, 0, 0)
        jg_history = zeros(Float64, 0, 0)
    end

    for n in eachindex(time_grid)
        t_now = time_grid[n]

        activate_nodes!(
            jq,
            jg,
            active,
            boundary_values.jq_bc,
            boundary_values.jg_bc,
            lattice.t_start,
            t_now,
        )

        if store_history
            jq_history[:, n] .= jq
            jg_history[:, n] .= jg
        end

        if n == length(time_grid)
            break
        end

        dt = time_grid[n + 1] - t_now

        if method == :euler
            components = build_splitting_component_grid(
                delta_ell = lattice.delta_ell,
                n_nodes = lattice.n_nodes,
                mu = t_to_mu(t_now),
                order = order,
                alpha_s_Z_value = alpha_s_Z_value,
                alpha_s_loops = alpha_s_loops,
            )

            rhs = build_rhs(jq = jq, jg = jg, active = active, components = components)

            for i in eachindex(active)
                if active[i]
                    jq[i] += dt * rhs.rhs_q[i]
                    jg[i] += dt * rhs.rhs_g[i]
                end
            end

        elseif method == :rk2
            components_1 = build_splitting_component_grid(
                delta_ell = lattice.delta_ell,
                n_nodes = lattice.n_nodes,
                mu = t_to_mu(t_now),
                order = order,
                alpha_s_Z_value = alpha_s_Z_value,
                alpha_s_loops = alpha_s_loops,
            )

            k1 = build_rhs(jq = jq, jg = jg, active = active, components = components_1)

            jq_mid = copy(jq)
            jg_mid = copy(jg)

            for i in eachindex(active)
                if active[i]
                    jq_mid[i] += 0.5 * dt * k1.rhs_q[i]
                    jg_mid[i] += 0.5 * dt * k1.rhs_g[i]
                end
            end

            components_2 = build_splitting_component_grid(
                delta_ell = lattice.delta_ell,
                n_nodes = lattice.n_nodes,
                mu = t_to_mu(t_now + 0.5 * dt),
                order = order,
                alpha_s_Z_value = alpha_s_Z_value,
                alpha_s_loops = alpha_s_loops,
            )

            k2 = build_rhs(jq = jq_mid, jg = jg_mid, active = active, components = components_2)

            for i in eachindex(active)
                if active[i]
                    jq[i] += dt * k2.rhs_q[i]
                    jg[i] += dt * k2.rhs_g[i]
                end
            end

        else
            error("Unknown method: $method. Expected :euler or :rk2.")
        end
    end

    return (
        jq = jq,
        jg = jg,
        active = active,
        jq_bc = boundary_values.jq_bc,
        jg_bc = boundary_values.jg_bc,
        jq_history = jq_history,
        jg_history = jg_history,
        lattice = lattice,
        time_grid = collect(time_grid),
    )
end
