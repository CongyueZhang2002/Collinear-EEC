using Printf

if !isdefined(@__MODULE__, :solve_jet_rg)
    include(joinpath(@__DIR__, "numerical_rg.jl"))
end

const DEFAULT_MU_VALUES = Float64[10.0, 20.0, 50.0, 100.0, 200.0, 500.0]

function jet_boundary_coefficients(nf::Int64)
    zeta2 = pi^2 / 6.0
    zeta3 = 1.2020569031595942854
    zeta4 = pi^4 / 90.0

    jq1 = -(37.0 / 3.0) * CF
    jg1 = -(898.0 / 75.0) * CA - (14.0 / 25.0) * nf

    jq2 = (
        (
            152.0 * zeta4 - 478.0 * zeta3 - 106.0 * zeta2 +
            3498505.0 / 5184.0
        ) * CF^2 +
        (
            -76.0 * zeta4 + 280.0 * zeta3 + (1063.0 / 15.0) * zeta2 -
            164883727.0 / 324000.0
        ) * CA * CF +
        ((9.0 / 5.0) * zeta2 + 703847.0 / 24000.0) * CF * nf
    )

    jg2 = (
        (
            76.0 * zeta4 - (1054.0 / 5.0) * zeta3 -
            (2159.0 / 75.0) * zeta2 + 133639871.0 / 810000.0
        ) * CA^2 +
        (
            (44.0 / 5.0) * zeta3 - (127.0 / 25.0) * zeta2 +
            68111303.0 / 1620000.0
        ) * CA * nf +
        (
            4.0 * zeta3 + (14.0 / 5.0) * zeta2 -
            1528667.0 / 108000.0
        ) * CF * nf +
        (-(8.0 / 15.0) * zeta2 + 2344.0 / 1125.0) * nf^2
    )

    return (jq1 = jq1, jq2 = jq2, jg1 = jg1, jg2 = jg2)
end

"""
    build_physical_boundary_func(; order, nf_scheme, alpha_s_order)

Build the physical boundary condition with zero-log perturbative constants
through `order`, multiplied by `exp(-2.3b)` for quarks and `exp(-3.8b)` for
gluons.
"""
function build_physical_boundary_func(;
    order::Int64 = 2,
    nf_scheme::Symbol = :VFNS,
    alpha_s_order::Int64 = 4,
)
    0 <= order <= 2 ||
        throw(ArgumentError("The physical boundary supports order = 0, 1, or 2."))

    function boundary_func(; b, bstar, mu_start)
        nf = nf_func(mu_start; scheme = nf_scheme)
        a_s = alpha_s_func_numerical(
            mu_f = mu_start,
            order = alpha_s_order,
            nf_scheme = nf_scheme,
        ) / (4.0 * pi)
        coefficient = jet_boundary_coefficients(nf)
        jq_pert = 1.0
        jg_pert = 1.0
        if order >= 1
            jq_pert += a_s * coefficient.jq1
            jg_pert += a_s * coefficient.jg1
        end
        if order >= 2
            jq_pert += a_s^2 * coefficient.jq2
            jg_pert += a_s^2 * coefficient.jg2
        end
        return (
            jq_pert * exp(-2.3 * b),
            jg_pert * exp(-3.8 * b),
        )
    end

    return boundary_func
end

"""
    build_evolved_interpolator(; n_nodes, b_min, b_max, bstar_func,
                               boundary_func, order, nf_scheme,
                               alpha_s_order, method, closure_check)

Build a callable quark-gluon RG solution. The grid is constructed internally
from `n_nodes`, `b_min`, `b_max`, and `bstar_func`. If `boundary_func` is
omitted, the physical perturbative-times-nonperturbative boundary is used.
"""
function build_evolved_interpolator(;
    n_nodes::Int64 = 1000,
    b_min::Float64 = 0.001,
    b_max::Float64 = 50.0,
    bstar_func::Function = b -> b / sqrt(1.0 + b^2 / b0^2),
    boundary_func::Union{Nothing,Function} = nothing,
    order::Int64 = 2,
    nf_scheme::Symbol = :VFNS,
    alpha_s_order::Int64 = 4,
    method::Symbol = :rk2,
    closure_check::Symbol = :error,
)
    if boundary_func === nothing
        boundary_func = build_physical_boundary_func(
            order = order,
            nf_scheme = nf_scheme,
            alpha_s_order = alpha_s_order,
        )
    end

    return solve_jet_rg(
        n_nodes = n_nodes,
        b_min = b_min,
        b_max = b_max,
        bstar_func = bstar_func,
        boundary_func = boundary_func,
        order = order,
        nf_scheme = nf_scheme,
        alpha_s_order = alpha_s_order,
        method = method,
        closure_check = closure_check,
    )
end

function _mu_file_label(mu::Float64)
    return isinteger(mu) ? string(Int64(mu)) : replace(string(mu), "." => "p")
end

"""
    write_evolved_csvs(solution; output_dir, mu_values, b_values)

Sample a callable RG solution and write one `b,Jq,Jg` CSV file per scale.
"""
function write_evolved_csvs(
    solution::StepwiseRGSolution;
    output_dir::AbstractString = joinpath(@__DIR__, "data"),
    mu_values::AbstractVector{Float64} = DEFAULT_MU_VALUES,
    b_values::AbstractVector{Float64} = collect(exp.(
        range(log(0.01), log(5.0), length = 500),
    )),
)
    mkpath(output_dir)
    output_paths = String[]

    for mu in mu_values
        output_path = joinpath(
            output_dir,
            "jet_rg_mu_$(_mu_file_label(mu))_GeV.csv",
        )
        open(output_path, "w") do io
            println(io, "b,Jq,Jg")
            for b in b_values
                value = solution(b, mu)
                @printf(io, "%.16e,%.16e,%.16e\n", b, value.jq, value.jg)
            end
        end
        push!(output_paths, output_path)
    end

    return output_paths
end

if abspath(PROGRAM_FILE) == @__FILE__
    n_nodes = isempty(ARGS) ? 1000 : parse(Int64, only(ARGS))
    timed_solution = @timed build_evolved_interpolator(n_nodes = n_nodes)
    output_paths = write_evolved_csvs(timed_solution.value)

    @printf(
        "Built %d x %d NNLO VFNS interpolator in %.6f seconds.\n",
        n_nodes,
        n_nodes,
        timed_solution.time,
    )
    for output_path in output_paths
        println("Wrote ", output_path)
    end
end
