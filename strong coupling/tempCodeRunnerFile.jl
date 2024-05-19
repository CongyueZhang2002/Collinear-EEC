# up to 4-loop beta function and strong coupling https://arxiv.org/pdf/2403.04077

using Distributed
using DifferentialEquations
include("constants.jl")

# μf: final scale
# μi: initial scale
# αs: alpha_s at initial scale
# order: 1, 2, 3, 4
# nf: active quarks

@everywhere function β0_func(nf) 
    value = (11*CA - 2*nf)/12
    return 4*value
end

@everywhere function β1_func(nf) 
    value = (17*CA^2 - 5*CA*nf - 3*CF*nf)/24
    return 4^2*value
end

@everywhere function β2_func(nf)
    value = (
    2857/54*CA^3 
    - (1415/54*CA^2 + 205/18*CA*CF - CF^2)*nf 
    + (79/54*CA + 11/9*CF)*nf^2
    )/64
    return 4^3*value
end

@everywhere function β3_func(nf)
    value = (
    1093/186624*nf^3
    + (809*z3/2592 + 50065/41472)*nf^2
    + (- 1627*z3/1728 - 1078361/41472)*nf
    + 891/64*z3 + 149753/1536
    )
    return 4^4*value
end

print(β2_func(5))