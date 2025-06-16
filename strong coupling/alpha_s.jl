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

    if nf == 4
        value = 8.333336
    elseif nf == 5
        value = 7.66667
    end

    return value 
end

@everywhere function β1_func(nf) 

    if nf == 4
        value = 51.3332
    elseif nf == 5
        value = 38.6665
    end

    return value
end

@everywhere function β2_func(nf)

    if nf == 4
        value = 406.352
    elseif nf == 5
        value = 180.908
    end

    return value
end

@everywhere function β3_func(nf)

    if nf == 4
        value = 8035.21
    elseif nf == 5
        value = 4826.18
    end

    return value
end

@everywhere function β4_func(nf)

    if nf == 4
        value = 58311.3
    elseif nf == 5
        value = 15471.7
    end

    return value
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

#usual resummed one
@everywhere function alpha_s_func(μf)

    order = nloops_αs

    β0 = β0_func(5)
    β1 = β1_func(5)
    β2 = β2_func(5)
    β3 = β3_func(5)

    μ_ini = 91.2
    l = 1 + β0*αs_Z*log(μf^2/μ_ini^2)/(4*π)
        
    order1 = 1 
    order2 = - αs_Z/(4π*l)*β1/β0*log(l)
    order3 = (αs_Z/(4π*l))^2 * (
            (β1/β0)^2 * (log(l)^2 - log(l) + l - 1) - β2/β0 * (l - 1)
            )
    order4 = (αs_Z/(4π*l))^3 * (
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

        αs_Z_final = αs_Z/l*total
    
    return αs_Z_final

end

@everywhere function alpha_s_freeze(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(5)
    β1 = β1_func(5)
    β2 = β2_func(5)
    β3 = β3_func(5)

    if μf > 1
        μfreeze = μf
    else
        μfreeze = 1
    end
    
    l = 1 + β0*αs*log(μfreeze^2/μi^2)/(4*π)
        
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

@everywhere function alpha_s_resum(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)

    l = 1 + β0*αs*log(μf^2/μi^2)/(4*π)
        
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

    return αs*total
end

function alpha_s_func_pert_freeze(; μf::Float64, μi::Float64, αs::Float64, order::Int64, nf::Int64)

    a0 = αs/π

    if μf > 1
        μfreeze = μf
    else
        μfreeze = 1
    end

    L = log(μfreeze/μi)

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
