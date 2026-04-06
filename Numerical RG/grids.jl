struct LatticeGrid
    n_nodes::Int64
    b_min::Float64
    b_max::Float64
    delta_ell::Float64
    ell_grid::Vector{Float64}
    b_grid::Vector{Float64}
    bstar_grid::Vector{Float64}
    mu_i_grid::Vector{Float64}
    t_grid::Vector{Float64}
end

function build_lattice_grid(;
    n_nodes::Int64,
    b_min::Float64,
    b_max::Float64,
    bstar_func::Function,
)

    if n_nodes < 2
        error("n_nodes must be at least 2.")
    end
    if b_min <= 0.0
        error("b_min must be positive.")
    end
    if b_max < b_min
        error("b_max must be greater than or equal to b_min.")
    end

    ell_min = log(b_min)
    ell_max = log(b_max)
    delta_ell = (ell_max - ell_min) / (n_nodes - 1)
    ell_grid = collect(range(ell_max, ell_min, length = n_nodes))
    b_grid = exp.(ell_grid)
    bstar_grid = bstar_func.(b_grid)
    mu_i_grid = b0 ./ bstar_grid
    t_grid = log.(mu_i_grid.^2)

    return LatticeGrid(
        n_nodes,
        b_min,
        b_max,
        delta_ell,
        ell_grid,
        b_grid,
        bstar_grid,
        mu_i_grid,
        t_grid,
    )
end

function build_boundary_values(; grid::LatticeGrid, boundary_func::Function)

    jq_bc = zeros(Float64, grid.n_nodes)
    jg_bc = zeros(Float64, grid.n_nodes)

    for i in 1:grid.n_nodes
        jq_bc[i], jg_bc[i] = boundary_func(
            b = grid.b_grid[i],
            bstar = grid.bstar_grid[i],
            mu_start = grid.mu_i_grid[i],
        )
    end

    return (jq_bc = jq_bc, jg_bc = jg_bc)
end
