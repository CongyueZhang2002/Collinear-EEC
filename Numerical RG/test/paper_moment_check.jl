using Printf
using Test

include(joinpath(@__DIR__, "..", "numerical_rg.jl"))
include(joinpath(@__DIR__, "paper_moment_tests.jl"))

function print_moment_errors(label, actual, expected)
    for channel in propertynames(expected)
        observed = getproperty(actual, channel)
        reference = getproperty(expected, channel)
        absolute_error = abs(observed - reference)
        relative_error = absolute_error / max(abs(reference), eps())
        @printf(
            "%-12s %-2s actual=% .12e reference=% .12e abs=%.3e rel=%.3e\n",
            label,
            String(channel),
            observed,
            reference,
            absolute_error,
            relative_error,
        )
    end
end

nf = nf_func(10.0; scheme = :VFNS)

for (index, expected) in enumerate(paper_N3_references(nf))
    loop_order = index - 1
    print_moment_errors(
        "N=3 L=$loop_order",
        splitting_moment_N3(loop_order; nf = nf),
        expected,
    )
end

logarithmic = paper_log_moment_references(nf)
print_moment_errors(
    "LO log(x)",
    splitting_moment_N3(0; log_power = 1, nf = nf),
    logarithmic.lo_first,
)
print_moment_errors(
    "LO log(x)^2",
    splitting_moment_N3(0; log_power = 2, nf = nf),
    logarithmic.lo_second,
)
print_moment_errors(
    "NLO log(x)",
    splitting_moment_N3(1; log_power = 1, nf = nf),
    logarithmic.nlo_first,
)
