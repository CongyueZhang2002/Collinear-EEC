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
        throw(ArgumentError("n_nodes must be at least 2."))
    end
    if !isfinite(b_min) || b_min <= 0.0
        throw(ArgumentError("b_min must be finite and positive."))
    end
    if !isfinite(b_max) || b_max <= b_min
        throw(ArgumentError("b_max must be finite and greater than b_min."))
    end

    ell_min = log(b_min)
    ell_max = log(b_max)
    delta_ell = (ell_max - ell_min) / (n_nodes - 1)
    if !isfinite(delta_ell) || delta_ell <= 0.0
        throw(ArgumentError("b_min and b_max must define a representable log-space interval."))
    end

    ell_grid = collect(range(ell_max, ell_min, length = n_nodes))
    b_grid = exp.(ell_grid)
    bstar_grid = Float64.(bstar_func.(b_grid))

    if !all(isfinite, bstar_grid) || !all(>(0.0), bstar_grid)
        throw(DomainError(bstar_grid, "bstar_func must return finite positive values."))
    end

    mu_i_grid = b0 ./ bstar_grid
    t_grid = log.(mu_i_grid.^2)

    if !all(isfinite, mu_i_grid) || !all(isfinite, t_grid)
        throw(DomainError(mu_i_grid, "The initial scales and RG times must be finite."))
    end
    if !all(>(0.0), diff(t_grid))
        throw(ArgumentError("bstar_func must produce a strictly increasing t_grid."))
    end

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

        if !isfinite(jq_bc[i]) || !isfinite(jg_bc[i])
            throw(DomainError(
                (jq_bc[i], jg_bc[i]),
                "boundary_func must return finite quark and gluon values.",
            ))
        end
    end

    return (jq_bc = jq_bc, jg_bc = jg_bc)
end
