using DifferentialEquations

if !isdefined(@__MODULE__, :β0_func)
    include("anomalous dims.jl")
end

const NF_SCHEMES = (:VFNS, :nf5)

function _check_nf_scheme(nf_scheme::Symbol)
    nf_scheme in NF_SCHEMES ||
        throw(ArgumentError("nf_scheme must be :VFNS or :nf5."))
    return nf_scheme
end

"""
    nf_func(mu; scheme=:VFNS)

Return the active flavor count at scale `mu`. `:VFNS` switches at the
temporary heavy-quark masses in `constants.jl`; `:nf5` always returns five.
"""
function nf_func(mu::Real; scheme::Symbol = :VFNS)
    isfinite(mu) && mu > 0 || throw(DomainError(mu, "mu must be finite and positive."))
    _check_nf_scheme(scheme)

    scheme == :nf5 && return 5
    mu >= mt && return 6
    mu >= mb && return 5
    mu >= mc && return 4
    return 3
end

function m_func(nf::Integer)
    nf == 4 && return mc
    nf == 5 && return mb
    nf == 6 && return mt
    throw(ArgumentError("m_func requires nf = 4, 5, or 6."))
end

function ΛQCD_func(nf::Integer)
    nf == 3 && return ΛQCD3
    nf == 4 && return ΛQCD4
    nf == 5 && return ΛQCD5
    throw(ArgumentError("ΛQCD_func currently has values only for nf = 3, 4, or 5."))
end

function beta_coefficient(loop_index::Integer, nf::Integer)
    loop_index == 0 && return β0_func(nf)
    loop_index == 1 && return β1_func(nf)
    loop_index == 2 && return β2_func(nf)
    loop_index == 3 && return β3_func(nf)
    loop_index == 4 && return β4_func(nf)
    throw(ArgumentError("loop_index must be between 0 and 4."))
end

"""
    beta_alpha_s(alpha_s; order, nf)

Return `d alpha_s / d log(mu)` through `order` beta-function loops. This is a
low-level coefficient function; scale-dependent callers obtain `nf` from
`nf_func`.
"""
function beta_alpha_s(alpha_s::Real; order::Integer, nf::Integer)
    1 <= order <= 5 || throw(ArgumentError("order must be between 1 and 5."))
    isfinite(alpha_s) && alpha_s > 0 ||
        throw(DomainError(alpha_s, "alpha_s must be finite and positive."))
    _check_nf(nf)

    a_s = Float64(alpha_s) / (4.0 * pi)
    series = 0.0
    a_s_power = a_s
    for loop_index in 0:(order - 1)
        series += beta_coefficient(loop_index, nf) * a_s_power
        a_s_power *= a_s
    end
    return -2.0 * Float64(alpha_s) * series
end

# Compatibility with the notation used by the copied TMD beta-function code.
function β_func(;
    αs::Float64,
    order::Int64,
    mu::Float64,
    nf_scheme::Symbol = :VFNS,
)
    nf = nf_func(mu; scheme = nf_scheme)
    return beta_alpha_s(αs; order = order, nf = nf)
end

function _validate_alpha_s_inputs(;
    mu_ref::Real,
    alpha_s_ref::Real,
    order::Integer,
    nf_scheme::Symbol,
    reltol::Real,
    abstol::Real,
)
    isfinite(mu_ref) && mu_ref > 0 ||
        throw(ArgumentError("mu_ref must be finite and positive."))
    isfinite(alpha_s_ref) && alpha_s_ref > 0 ||
        throw(ArgumentError("alpha_s_ref must be finite and positive."))
    1 <= order <= 5 || throw(ArgumentError("order must be between 1 and 5."))
    _check_nf_scheme(nf_scheme)
    isfinite(reltol) && reltol > 0 ||
        throw(ArgumentError("reltol must be finite and positive."))
    isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("abstol must be finite and positive."))
    return nothing
end

function _threshold_tstops(
    t_ref::Float64,
    t_target::Float64,
    nf_scheme::Symbol,
)
    nf_scheme == :nf5 && return Float64[]

    t_low, t_high = minmax(t_ref, t_target)
    stops = Float64[
        2.0 * log(mass) for mass in (mc, mb, mt)
        if t_low < 2.0 * log(mass) < t_high
    ]
    sort!(stops)
    t_target < t_ref && reverse!(stops)
    return stops
end

function _solve_alpha_s_segment(
    t_ref::Float64,
    t_target::Float64,
    alpha_s_ref::Float64,
    order::Int64,
    nf_scheme::Symbol,
    reltol::Float64,
    abstol::Float64,
)
    function rhs(alpha_s, _, t)
        mu = exp(0.5 * t)
        nf = nf_func(mu; scheme = nf_scheme)
        return 0.5 * beta_alpha_s(alpha_s; order = order, nf = nf)
    end

    problem = ODEProblem(rhs, alpha_s_ref, (t_ref, t_target))
    solution = solve(
        problem,
        Tsit5();
        reltol = reltol,
        abstol = abstol,
        tstops = _threshold_tstops(t_ref, t_target, nf_scheme),
        verbose = false,
    )
    endpoint_tolerance = 16.0 * eps(max(abs(t_ref), abs(t_target), 1.0))
    if !isapprox(solution.t[end], t_target; rtol = 0.0, atol = endpoint_tolerance)
        throw(DomainError(
            t_target,
            "Numerical alpha_s evolution did not reach the requested scale.",
        ))
    end
    return solution
end

"""
    alpha_s_grid(t_grid; mu_ref=MZ, alpha_s_ref=DEFAULT_ALPHA_S_MZ,
                 order, nf_scheme=:VFNS)

Numerically evolve the coupling from `(mu_ref, alpha_s_ref)` and evaluate it
at `t_grid = log(mu^2)`. In `:VFNS`, the beta function changes at `mc`, `mb`,
and `mt`, with a continuous coupling at each threshold. `:nf5` uses five
active flavors at every scale.
"""
function alpha_s_grid(
    t_grid::AbstractVector{<:Real};
    mu_ref::Real = MZ,
    alpha_s_ref::Real = DEFAULT_ALPHA_S_MZ,
    order::Integer,
    nf_scheme::Symbol = :VFNS,
    reltol::Real = 1.0e-10,
    abstol::Real = 1.0e-12,
)
    isempty(t_grid) && throw(ArgumentError("t_grid must not be empty."))
    all(isfinite, t_grid) || throw(ArgumentError("t_grid must be finite."))
    _validate_alpha_s_inputs(
        mu_ref = mu_ref,
        alpha_s_ref = alpha_s_ref,
        order = order,
        nf_scheme = nf_scheme,
        reltol = reltol,
        abstol = abstol,
    )

    t_values = Float64.(t_grid)
    t_ref = 2.0 * log(Float64(mu_ref))
    t_min, t_max = extrema(t_values)
    lower_solution = t_min < t_ref ? _solve_alpha_s_segment(
        t_ref,
        t_min,
        Float64(alpha_s_ref),
        Int64(order),
        nf_scheme,
        Float64(reltol),
        Float64(abstol),
    ) : nothing
    upper_solution = t_max > t_ref ? _solve_alpha_s_segment(
        t_ref,
        t_max,
        Float64(alpha_s_ref),
        Int64(order),
        nf_scheme,
        Float64(reltol),
        Float64(abstol),
    ) : nothing

    values = similar(t_values)
    for i in eachindex(t_values)
        t = t_values[i]
        values[i] = if t < t_ref
            lower_solution(t)
        elseif t > t_ref
            upper_solution(t)
        else
            Float64(alpha_s_ref)
        end
    end

    if !all(isfinite, values) || !all(>(0.0), values)
        throw(DomainError(values, "Numerical alpha_s evolution became non-finite or non-positive."))
    end
    return values
end

function alpha_s_func_numerical(;
    mu_f::Real,
    mu_ref::Real = MZ,
    alpha_s_ref::Real = DEFAULT_ALPHA_S_MZ,
    order::Integer,
    nf_scheme::Symbol = :VFNS,
    reltol::Real = 1.0e-10,
    abstol::Real = 1.0e-12,
)
    isfinite(mu_f) && mu_f > 0 ||
        throw(ArgumentError("mu_f must be finite and positive."))
    return only(alpha_s_grid(
        [2.0 * log(Float64(mu_f))];
        mu_ref = mu_ref,
        alpha_s_ref = alpha_s_ref,
        order = order,
        nf_scheme = nf_scheme,
        reltol = reltol,
        abstol = abstol,
    ))
end

function αs_func(
    mu::Real;
    mu_ref::Real = MZ,
    alpha_s_ref::Real = DEFAULT_ALPHA_S_MZ,
    order::Integer,
    nf_scheme::Symbol = :VFNS,
    reltol::Real = 1.0e-10,
    abstol::Real = 1.0e-12,
)
    return alpha_s_func_numerical(
        mu_f = mu,
        mu_ref = mu_ref,
        alpha_s_ref = alpha_s_ref,
        order = order,
        nf_scheme = nf_scheme,
        reltol = reltol,
        abstol = abstol,
    )
end

function alpha_qed(Q::Real)
    Q > 0 || throw(DomainError(Q, "Q must be positive."))

    mtau = 1.777
    alpha_em_mc = 0.007470670604
    alpha_em_mtau = 0.007491983578
    alpha_em_mb = 0.007557200440
    alpha_em_mz = 0.007815247548
    mc_qed = 1.2700
    mtau_qed = 1.777
    mb_qed = 4.180
    mz_qed = MZ

    eq23 = 6.0 / 9.0
    eq24 = 10.0 / 9.0
    eq25 = 11.0 / 9.0

    if Q <= mc_qed
        b0_qed = (2 + 3 * eq23) / (3 * pi)
        return alpha_em_mc / (1 - b0_qed * alpha_em_mc * log(Q^2 / mc_qed^2))
    elseif Q <= mtau
        b0_qed = (2 + 3 * eq24) / (3 * pi)
        return alpha_em_mc / (1 - b0_qed * alpha_em_mc * log(Q^2 / mc_qed^2))
    elseif Q <= mb_qed
        b0_qed = (3 + 3 * eq24) / (3 * pi)
        return alpha_em_mtau / (
            1 - b0_qed * alpha_em_mtau * log(Q^2 / mtau_qed^2)
        )
    elseif Q <= mz_qed
        b0_qed = (3 + 3 * eq25) / (3 * pi)
        return alpha_em_mb / (1 - b0_qed * alpha_em_mb * log(Q^2 / mb_qed^2))
    end

    b0_qed = (3 + 3 * eq25) / (3 * pi)
    return alpha_em_mz / (1 - b0_qed * alpha_em_mz * log(Q^2 / mz_qed^2))
end
