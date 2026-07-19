using Printf
using Test

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))
include(joinpath(@__DIR__, "paper_moment_tests.jl"))

function moment_errors(nf, label, actual, expected)
    return [
        begin
            observed = getproperty(actual, channel)
            reference = getproperty(expected, channel)
            absolute = abs(observed - reference)
            relative = absolute / max(abs(reference), eps(Float64))
            (
                nf = nf,
                label = label,
                channel = channel,
                observed = observed,
                reference = reference,
                absolute = absolute,
                relative = relative,
            )
        end
        for channel in propertynames(expected)
    ]
end

function maximum_error(records, field)
    return reduce(records) do largest, record
        getproperty(record, field) > getproperty(largest, field) ? record : largest
    end
end

all_records = NamedTuple[]

println("Maximum absolute differences from arXiv:1905.01310")
println("nf   LO N=3       NLO N=3      NNLO N=3     logarithmic")

for nf in 3:6
    references = paper_N3_references(nf)
    nf_records = NamedTuple[]
    n3_maxima = Float64[]

    for loop_order in 0:2
        records = moment_errors(
            nf,
            "N=3 L=$loop_order",
            splitting_moment_N3(loop_order; nf = nf),
            references[loop_order + 1],
        )
        append!(nf_records, records)
        push!(n3_maxima, maximum(record.absolute for record in records))
    end

    logarithmic = paper_log_moment_references(nf)
    log_records = vcat(
        moment_errors(
            nf,
            "LO log(x)",
            splitting_moment_N3(0; log_power = 1, nf = nf),
            logarithmic.lo_first,
        ),
        moment_errors(
            nf,
            "LO log(x)^2",
            splitting_moment_N3(0; log_power = 2, nf = nf),
            logarithmic.lo_second,
        ),
        moment_errors(
            nf,
            "NLO log(x)",
            splitting_moment_N3(1; log_power = 1, nf = nf),
            logarithmic.nlo_first,
        ),
    )
    append!(nf_records, log_records)
    append!(all_records, nf_records)

    @printf(
        "%d    %.6e   %.6e   %.6e   %.6e\n",
        nf,
        n3_maxima[1],
        n3_maxima[2],
        n3_maxima[3],
        maximum(record.absolute for record in log_records),
    )
end

worst_absolute = maximum_error(all_records, :absolute)
worst_relative = maximum_error(all_records, :relative)

println()
@printf(
    "Worst absolute: nf=%d %s %s observed=% .12e reference=% .12e abs=%.12e rel=%.12e\n",
    worst_absolute.nf,
    worst_absolute.label,
    String(worst_absolute.channel),
    worst_absolute.observed,
    worst_absolute.reference,
    worst_absolute.absolute,
    worst_absolute.relative,
)
@printf(
    "Worst relative: nf=%d %s %s observed=% .12e reference=% .12e abs=%.12e rel=%.12e\n",
    worst_relative.nf,
    worst_relative.label,
    String(worst_relative.channel),
    worst_relative.observed,
    worst_relative.reference,
    worst_relative.absolute,
    worst_relative.relative,
)
