struct LatticeGrid
    n_nodes::Int64
    ell_max::Float64
    delta_ell::Float64
    bmax::Float64
    ell::Vector{Float64}
    b::Vector{Float64}
    bstar::Vector{Float64}
    mu_start::Vector{Float64}
    t_start::Vector{Float64}
end

function build_lattice_grid(; n_nodes::Int64, ell_max::Float64, delta_ell::Float64, bmax::Float64)

    ell = [ell_max - (i - 1) * delta_ell for i in 1:n_nodes]
    b = exp.(ell)
    bstar = b ./ sqrt.(1 .+ (b ./ bmax).^2)
    mu_start = b0 ./ bstar
    t_start = log.(mu_start.^2)

    return LatticeGrid(
        n_nodes,
        ell_max,
        delta_ell,
        bmax,
        ell,
        b,
        bstar,
        mu_start,
        t_start,
    )
end

function build_time_grid(; t_start::AbstractVector{Float64}, delta_t::Float64, t_final::Float64)

    t_min = minimum(t_start)

    if t_final < t_min
        error("t_final must be at least the smallest activation time.")
    end

    regular_grid = collect(t_min:delta_t:t_final)

    if isempty(regular_grid) || regular_grid[end] < t_final
        push!(regular_grid, t_final)
    end

    merged_grid = sort!(unique(vcat(regular_grid, collect(t_start), [t_final])))

    return merged_grid
end

function build_boundary_values(; grid::LatticeGrid, boundary_func::Function)

    jq_bc = zeros(Float64, grid.n_nodes)
    jg_bc = zeros(Float64, grid.n_nodes)

    for i in 1:grid.n_nodes
        jq_bc[i], jg_bc[i] = boundary_func(
            b = grid.b[i],
            bstar = grid.bstar[i],
            mu_start = grid.mu_start[i],
            t_start = grid.t_start[i],
        )
    end

    return (jq_bc = jq_bc, jg_bc = jg_bc)
end
