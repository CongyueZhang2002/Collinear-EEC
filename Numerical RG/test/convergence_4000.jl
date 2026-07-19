using Printf

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))

bstar_func(b) = b / sqrt(1.0 + b^2 / b0^2)

boundary_func(; b, bstar, mu_start) = (
    exp(-2.0 * b),
    0.5 * exp(-3.0 * b),
)

function build_timed_solution(n_nodes::Int64)
    GC.gc()
    result = @timed solve_jet_rg(
        n_nodes = n_nodes,
        b_min = 0.001,
        b_max = 30.0,
        bstar_func = bstar_func,
        boundary_func = boundary_func,
        order = 2,
        method = :rk2,
    )
    @printf("SOLVE n_nodes=%d seconds=%.6f\n", n_nodes, result.time)
    flush(stdout)
    return result.value
end

function sample_solution(
    solution::StepwiseRGSolution,
    ell_validation::Vector{Float64},
    t_validation::Vector{Float64},
)

    jq = Matrix{Float64}(undef, length(ell_validation), length(t_validation))
    jg = similar(jq)

    Threads.@threads :static for n in eachindex(t_validation)
        mu = exp(0.5 * t_validation[n])
        for i in eachindex(ell_validation)
            b = clamp(
                exp(ell_validation[i]),
                solution.lattice.b_min,
                solution.lattice.b_max,
            )
            value = solution(b, mu)
            jq[i, n] = value.jq
            jg[i, n] = value.jg
        end
    end

    return jq, jg
end

function error_metrics(values::AbstractArray{Float64}, reference::AbstractArray{Float64})
    difference = values .- reference
    reference_scale = maximum(abs, reference)
    relative_floor = 1e-6 * reference_scale
    relative_errors = Float64[]

    for i in eachindex(values, reference)
        if abs(reference[i]) >= relative_floor
            push!(relative_errors, abs(difference[i] / reference[i]))
        end
    end

    sort!(relative_errors)
    p95_index = max(1, ceil(Int64, 0.95 * length(relative_errors)))

    return (
        normalized_rms = sqrt(sum(abs2, difference) / sum(abs2, reference)),
        normalized_max = maximum(abs, difference) / reference_scale,
        relative_p95 = relative_errors[p95_index],
    )
end

function print_metrics(label::String, metrics)
    @printf(
        "ERROR %s normalized_rms=%.8e normalized_max=%.8e relative_p95=%.8e\n",
        label,
        metrics.normalized_rms,
        metrics.normalized_max,
        metrics.relative_p95,
    )
end

# Compile all solver and interpolation paths before timing.
warmup_solution = solve_jet_rg(
    n_nodes = 12,
    b_min = 0.1,
    b_max = 1.0,
    bstar_func = bstar_func,
    boundary_func = boundary_func,
    order = 2,
    method = :rk2,
)
warmup_solution(0.5, sqrt(
    warmup_solution.lattice.mu_i_grid[1] *
    warmup_solution.lattice.mu_i_grid[end]
))
warmup_solution = nothing

println("Julia threads: ", Threads.nthreads())
solution_4000 = build_timed_solution(4000)
solution_2000 = build_timed_solution(2000)
solution_1000 = build_timed_solution(1000)

n_validation = 301
ell_validation = collect(range(log(30.0), log(0.001), length = n_validation))
t_validation = collect(range(
    solution_4000.time_grid[1],
    solution_4000.time_grid[end],
    length = n_validation,
))

jq_4000, jg_4000 = sample_solution(solution_4000, ell_validation, t_validation)
jq_2000, jg_2000 = sample_solution(solution_2000, ell_validation, t_validation)
jq_1000, jg_1000 = sample_solution(solution_1000, ell_validation, t_validation)

jq_2000_metrics = error_metrics(jq_2000, jq_4000)
jg_2000_metrics = error_metrics(jg_2000, jg_4000)
jq_1000_metrics = error_metrics(jq_1000, jq_2000)
jg_1000_metrics = error_metrics(jg_1000, jg_2000)

print_metrics("full_jq_2000_vs_4000", jq_2000_metrics)
print_metrics("full_jg_2000_vs_4000", jg_2000_metrics)
print_metrics("full_jq_1000_vs_2000", jq_1000_metrics)
print_metrics("full_jg_1000_vs_2000", jg_1000_metrics)

jq_order = log(
    jq_1000_metrics.normalized_rms / jq_2000_metrics.normalized_rms,
) / log(2.0)
jg_order = log(
    jg_1000_metrics.normalized_rms / jg_2000_metrics.normalized_rms,
) / log(2.0)
@printf("ORDER jq=%.6f jg=%.6f\n", jq_order, jg_order)

csv_path = joinpath(@__DIR__, "convergence_4000_vs_2000_by_mu.csv")
selected_indices = round.(Int64, range(1, n_validation, length = 11))

open(csv_path, "w") do csv
    println(
        csv,
        "mu,jq_rms_percent,jg_rms_percent,jq_max_percent,jg_max_percent",
    )

    for n in eachindex(t_validation)
        mu = exp(0.5 * t_validation[n])
        jq_metrics = error_metrics(@view(jq_2000[:, n]), @view(jq_4000[:, n]))
        jg_metrics = error_metrics(@view(jg_2000[:, n]), @view(jg_4000[:, n]))

        @printf(
            csv,
            "%.12g,%.10e,%.10e,%.10e,%.10e\n",
            mu,
            100.0 * jq_metrics.normalized_rms,
            100.0 * jg_metrics.normalized_rms,
            100.0 * jq_metrics.normalized_max,
            100.0 * jg_metrics.normalized_max,
        )

        if n in selected_indices
            @printf(
                "BY_MU mu=%.8g jq_rms=%.6f%% jg_rms=%.6f%%\n",
                mu,
                100.0 * jq_metrics.normalized_rms,
                100.0 * jg_metrics.normalized_rms,
            )
        end
    end
end

println("WROTE ", csv_path)
