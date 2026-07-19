struct StepwiseRGSolution
    jq::Vector{Float64}
    jg::Vector{Float64}
    active::BitVector
    jq_bc::Vector{Float64}
    jg_bc::Vector{Float64}
    jq_history::Matrix{Float64}
    jg_history::Matrix{Float64}
    lattice::LatticeGrid
    time_grid::Vector{Float64}
end

function _ell_cell(lattice::LatticeGrid, b::Real)
    b_value = Float64(b)
    if !isfinite(b_value) || !(lattice.b_min <= b_value <= lattice.b_max)
        throw(DomainError(
            b,
            "b must lie in [$(lattice.b_min), $(lattice.b_max)].",
        ))
    end

    ell = log(b_value)
    ell_position = (lattice.ell_grid[1] - ell) / lattice.delta_ell
    ell_position = clamp(ell_position, 0.0, lattice.n_nodes - 1.0)
    b_index = min(floor(Int64, ell_position) + 1, lattice.n_nodes - 1)
    b_fraction = ell_position - (b_index - 1)

    return b_index, b_fraction
end

"""
    mu_start(solution, b)

Return the piecewise-log-linear representation of the boundary scale `b0 / bstar(b)`.
"""
function mu_start(solution::StepwiseRGSolution, b::Real)
    b_index, b_fraction = _ell_cell(solution.lattice, b)
    t_start = (
        (1.0 - b_fraction) * solution.lattice.t_grid[b_index] +
        b_fraction * solution.lattice.t_grid[b_index + 1]
    )
    return exp(0.5 * t_start)
end

function _interpolation_cell(solution::StepwiseRGSolution, b::Real, mu::Real)
    n_nodes = solution.lattice.n_nodes
    expected_history_size = (n_nodes, n_nodes)

    if size(solution.jq_history) != expected_history_size ||
       size(solution.jg_history) != expected_history_size
        throw(ArgumentError(
            "Interpolation requires solve_stepwise_rg(..., store_history = true).",
        ))
    end

    mu_value = Float64(mu)
    if !isfinite(mu_value) || mu_value <= 0.0
        throw(DomainError(mu, "mu must be finite and positive."))
    end

    b_index, b_fraction = _ell_cell(solution.lattice, b)
    t = 2.0 * log(mu_value)
    t_min = solution.time_grid[1]
    t_final = solution.time_grid[end]
    time_tol = 64.0 * eps(Float64) * max(1.0, abs(t), abs(t_min), abs(t_final))

    if t < t_min - time_tol
        throw(DomainError(
            mu,
            "mu must be at least mu_min = $(exp(0.5 * t_min)).",
        ))
    end
    if t > t_final + time_tol
        throw(DomainError(
            mu,
            "mu must not exceed mu_max = $(exp(0.5 * t_final)).",
        ))
    end

    t = clamp(t, t_min, t_final)
    time_index = searchsortedlast(solution.time_grid, t)
    if time_index >= n_nodes
        time_index = n_nodes - 1
        time_fraction = 1.0
    else
        delta_t = solution.time_grid[time_index + 1] - solution.time_grid[time_index]
        time_fraction = (t - solution.time_grid[time_index]) / delta_t
        time_fraction = clamp(time_fraction, 0.0, 1.0)
    end

    on_boundary_cell = time_index == b_index
    above_boundary = false
    if on_boundary_cell
        fraction_tol = 128.0 * eps(Float64)
        if abs(time_fraction - b_fraction) <= fraction_tol
            time_fraction = b_fraction
        end
        above_boundary = time_fraction >= b_fraction
    end

    return (
        b_index = b_index,
        time_index = time_index,
        b_fraction = b_fraction,
        time_fraction = time_fraction,
        on_boundary_cell = on_boundary_cell,
        above_boundary = above_boundary,
    )
end

@inline function _interpolate_history(history::Matrix{Float64}, cell)
    i = cell.b_index
    n = cell.time_index
    b_fraction = cell.b_fraction
    time_fraction = cell.time_fraction

    if cell.on_boundary_cell
        if cell.above_boundary
            return (
                (1.0 - time_fraction) * history[i, n] +
                (time_fraction - b_fraction) * history[i, n + 1] +
                b_fraction * history[i + 1, n + 1]
            )
        end

        return (
            (1.0 - b_fraction) * history[i, n] +
            (b_fraction - time_fraction) * history[i + 1, n] +
            time_fraction * history[i + 1, n + 1]
        )
    end

    return (
        (1.0 - b_fraction) * (1.0 - time_fraction) * history[i, n] +
        b_fraction * (1.0 - time_fraction) * history[i + 1, n] +
        (1.0 - b_fraction) * time_fraction * history[i, n + 1] +
        b_fraction * time_fraction * history[i + 1, n + 1]
    )
end

function interpolate_jet(solution::StepwiseRGSolution, b::Real, mu::Real)
    cell = _interpolation_cell(solution, b, mu)
    return (
        jq = _interpolate_history(solution.jq_history, cell),
        jg = _interpolate_history(solution.jg_history, cell),
    )
end

"""
    solution(b, mu; component = :both)

Interpolate the evolved jet functions. The default result is `(jq = ..., jg = ...)`;
use `component = :q` or `:g` to return one scalar component. The valid scale
range is `solution.time_grid[1] <= log(mu^2) <= solution.time_grid[end]`.
"""
function (solution::StepwiseRGSolution)(
    b::Real,
    mu::Real;
    component::Symbol = :both,
)

    if !(component in (:both, :q, :g))
        throw(ArgumentError("component must be :both, :q, or :g."))
    end

    value = interpolate_jet(solution, b, mu)
    if component == :q
        return value.jq
    elseif component == :g
        return value.jg
    end
    return value
end
