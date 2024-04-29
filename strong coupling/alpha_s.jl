# up to 4-loop beta function and strong coupling https://arxiv.org/pdf/2403.04077

include("constants.jl")

# μf: final scale
# μi: initial scale
# αs: alpha_s at initial scale
# order: 1, 2, 3, 4
# nf: active quarks

function β0_func(nf) 
    value = (11*CA - 2*nf)/12
    return value
end

function β1_func(nf) 
    value = (17*CA^2 - 5*CA*nf - 3*CF*nf)/24
    return value
end

function β2_func(nf)
    value = (
    2857/54*CA^3 
    - (1415/54*CA^2 + 205/18*CA*CF - CF^2)*nf 
    + (79/54*CA + 11/9*CF)*nf^2
    )/64
    return value
end

function β3_func(nf)
    value = (
    1093/186624*nf^3
    + (809*z3/2592 + 50065/41472)*nf^2
    + (- 1627*z3/1728 - 1078361/41472)*nf
    + 891/64*z3 + 149753/1536
    )
    return value
end

function alpha_s_func(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)

    l = 1 + β0*αs*log(μf^2/μi^2)/π

    order1 = 1 
    order2 = - αs/(π*l)*β1/β0*log(l)
    order3 = (αs/(π*l))^2 * (
        (β1/β0)^2 * (log(l)^2 - log(l) + l - 1) - β2/β0 * (l - 1)
    )
    order4 = (αs/(π*l))^3 * (
        (β1/β0)^3 * (-log(l)^3 + 5/2*log(l)^2 - 2*(l - 1)*log(l) - (l-1)^2/2)
        + β1*β2/β0^2 * ((l - 1)*l + (2*l - 3)*log(l))
        + β3/β0 * (1 - l^2)/2
    )

    if order == 1 
        total = order1
    end
    if order == 2 
        total = order1 + order2
    end
    if order == 3 
        total = order1 + order2 + order3
    end
    if order == 4 
        total = order1 + order2 + order3 + order4
    end

    αs_final = αs/l*total

    return αs_final

end