# strong coupling https://arxiv.org/pdf/2403.04077
# 5-loop beta function https://arxiv.org/pdf/1701.01404

using Distributed
using DifferentialEquations
include("constants.jl")

# μf: final scale
# μi: initial scale
# αs: alpha_s at initial scale
# order(loops): 1, 2, 3, 4, 5
# nf: active quarks

@everywhere function β0_func(nf) 
    #value = (
    #  11 
    #- 0.666666*nf
    #)
    #return value
    return 7.66667
end

@everywhere function β1_func(nf) 
    #value = (
    #  102 
    #- 12.6667*nf
    #)
    #return value
    return 38.6665
end

@everywhere function β2_func(nf)
    #value = (
    #  1428.5 
    #- 279.611*nf 
    #+ 6.01852*nf^2
    #)
    #return value
    return 180.908
end

@everywhere function β3_func(nf)
    #value = (
    #  29243.0 
    #- 6946.29*nf
    #+ 405.089*nf^2
    #+ 1.49931*nf^3
    #)
    #return value
    return 4826.18
end

@everywhere function β4_func(nf)
    #value = (
    #  537148
    #- 186162*nf
    #+ 17567.8*nf^2
    #− 231.278*nf^3
    #− 1.84247*nf^4
    #)
    #return value
    return 15471.7
end

function β_func(; αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)
    β4 = β4_func(nf)

    order1 = αs/(4π)*β0
    order2 = (αs/(4π))^2*β1
    order3 = (αs/(4π))^3*β2
    order4 = (αs/(4π))^4*β3
    order5 = (αs/(4π))^4*β4

    if order == 1 
        total = order1
    elseif order == 2 
        total = order1 + order2
    elseif order == 3 
        total = order1 + order2 + order3
    elseif order == 4 
        total = order1 + order2 + order3 + order4
    elseif order == 5 
        total = order1 + order2 + order3 + order4 + order5
    end

    return -2*αs*total

end

@everywhere function alpha_s_func(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)

    l = 1 + β0*αs*log(μf^2/μi^2)/(4*π)
    if l < 10^(-3)
        l = 10^(-3)
    end
    
    order1 = 1 
    order2 = - αs/(4π*l)*β1/β0*log(l)
    order3 = (αs/(4π*l))^2 * (
        (β1/β0)^2 * (log(l)^2 - log(l) + l - 1) - β2/β0 * (l - 1)
    )
    order4 = (αs/(4π*l))^3 * (
        (β1/β0)^3 * (-log(l)^3 + 5/2*log(l)^2 - 2*(l - 1)*log(l) - (l-1)^2/2)
        + β1*β2/β0^2 * ((l - 1)*l + (2*l - 3)*log(l))
        + β3/β0 * (1 - l^2)/2
    )
 
    if order == 1 
        total = order1
    elseif order == 2 
        total = order1 + order2
    elseif order == 3 
        total = order1 + order2 + order3
    elseif order == 4 
        total = order1 + order2 + order3 + order4
    end

    αs_final = αs/l*total

    return αs_final

end

function alpha_s_func_schwartz(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)

    L = log(μf/μi)
    order1 = αs
    order2 = - αs^2/(2π)*β0*L
    order3 = αs^3/(8*π^2)*(-β1*L+2*β0^2*L^2)
    order4 = αs^4/(32*π^3)*(-β2*L+5*β0*β1*L^2-4*β0^3*L^3)
 
    if order == 1 
        total = order1 + order2
    end
    if order == 2 
        total = order1 + order2 + order3
    end
    if order == 3 
        total = order1 + order2 + order3 + order4
    end

    return total

end

function alpha_s_func_pert(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    a0 = αs/π
    L = log(μf/μi)

    β1 = -2*β0_func(nf)/4
    β2 = -2*β1_func(nf)/4^2
    β3 = -2*β2_func(nf)/4^3
    β4 = -2*β3_func(nf)/4^4
    β5 = -2*β4_func(nf)/4^5

    order1 = 1 + β1*L*a0 
    order2 = (β1^2*L^2 + β2*L)*a0^2
    order3 = (β1^3*L^3 + 5/2*β1*β2*L^2 + β3*L)*a0^3
    order4 = (β1^4*L^4 + 13/3*β1^2*β2*L^3 + (3/2*β2^2 + 3*β1*β3)*L^2 + β4*L)*a0^4
    order5 = (β1^5*L^5 + 77/12*β1^3*β2*L^4 + 1/6*(35*β1*β2^2 + 36*β1^2*β3)*L^3 + 7/2*(β1*β4 + β2*β3)*L^2 + β5*L)*a0^5

    if order == 0
        total = order1
    elseif order == 1 
        total = order1 + order2
    elseif order == 2 
        total = order1 + order2 + order3
    elseif order == 3 
        total = order1 + order2 + order3 + order4
    elseif order == 4 
        total = order1 + order2 + order3 + order4 + order5      
    end

    return αs*total
end

function alpha_s_func_numerical(; μf::Float64, μi::Float64, αs::Float64, order::Int64)

    function ode_function(α, p, μ)
        return 1/μ * β_func(αs=α[1], order=p, nf=5)
    end

    α0 = αs
    tspan = (μi, μf)  

    prob = ODEProblem(ode_function, α0, tspan, order)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    return sol[end][1] 
end

#function alpha_s_func_NNLL(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)
#
#    β0 = β0_func(nf)
#    β1 = β1_func(nf)
#    β2 = β2_func(nf)
#    β3 = β3_func(nf)
#
#    order1 = αs/(1+αs/(4π)*β0*log((μf/μi)^2))
#    order2 = -(order1)^2*β1/(4*π*β0)*log(αs/order1)
#    order3 = (order1)^3*(
#        β1^2/(4π*β0)^2*log(αs/order1)*(log(αs/order1)-1)
#        - (β1^2/(4π*β0)^2 - β2/(16*π^2*β0))*(1-log(αs/order1))
#    )
#
#    if order == 1 
#        total = order1
#    end
#    if order == 2 
#        total = order1 + order2
#    end
#    if order == 3 
#        total = order1 + order2 + order3
#    end
#    return total
#    
#end

#@everywhere function β0_func(nf) 
#    value = (11*CA - 2*nf)/12
#    return 4*value
#end

#@everywhere function β1_func(nf) 
#    value = (17*CA^2 - 5*CA*nf - 3*CF*nf)/24
#    return 4^2*value
#end

#@everywhere function β2_func(nf)
#    value = (
#    2857/54*CA^3 
#    - (1415/54*CA^2 + 205/18*CA*CF - CF^2)*nf 
#    + (79/54*CA + 11/9*CF)*nf^2
#    )/64
#    return 4^3*value
#end

#@everywhere function β3_func(nf)
#    value = (
#    1093/186624*nf^3
#    + (809*z3/2592 + 50065/41472)*nf^2
#    + (- 1627*z3/1728 - 1078361/41472)*nf
#    + 891/64*z3 + 149753/1536
#    )
#    return 4^4*value
#end