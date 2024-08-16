# Renormalon https://arxiv.org/pdf/2305.19311 (currently nf=5)
using HCubature
include("strong coupling\\alpha_s.jl")
include("strong coupling\\constants.jl")
include("fixed order/perturbation.jl")

function d1_func(μ,R,nf) 
    value = -8.357
    return value
end

function d2_func(μ,R,nf) 
    β0 = β0_func(5) 
    value = -72.443 - 16.713*β0*log(μ/R)
    return value
end

function γR0_func(nf) 
    return -8.357
end

function γR1_func(nf) 
    return 55.693
end

function ΣR_func(; μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    R1 = μ_ren
    μ = μ_ren

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    d1 = d1_func(μ,R1,nf) 
    d2 = d2_func(μ,R1,nf) 

    order1 = αs/(4π) * R1*d1
    order2 = (αs/(4π))^2 * R1*d2

    if order == 1 
        total = order1
    end
    if order == 2 
        total = order1 + order2
    end

    return total

end

function δR0_func(; μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    μ = μ_ren

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    d1 = d1_func(μ,R0,nf) 
    d2 = d2_func(μ,R0,nf) 

    order1 = αs/(4π) * R0*d1
    order2 = (αs/(4π))^2 * R0*d2

    if order == 1 
        total = order1
    end
    if order == 2 
        total = order1 + order2
    end

    return total

end

function Ω_Integral_func(; μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    R1 = μ_ren
    μ = μ_ren

    γR0 = γR0_func(nf)
    γR1 = γR1_func(nf)

    R_max = [R1]
    R_min = [R0]

    f1(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)/(4π)
    integral1, error1 = hcubature(f1, R_min, R_max)

    f2(x) = (alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)/(4π))^2
    integral2, error2 = hcubature(f2, R_min, R_max)

    order1 = γR0*integral1
    order2 = γR1*integral2

    if order == 1 
        total = order1
    end
    if order == 2 
        total = order1 + order2
    end

    return total

end

function renormalon_func(; z::Float64, Q::Float64, μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64, Ω1::Float64)

    ΣR = ΣR_func(μ_ren=μ_ren, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
    δR0 = δR0_func(μ_ren=μ_ren, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
    Ω_Integral = Ω_Integral_func(μ_ren=μ_ren, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    #total = (ΣR + Ω1 - δR0 - Ω_Integral)/(2*Q*(z*(1-z))^1.5)
    total = (ΣR- Ω_Integral)/(2*Q*(z*(1-z))^1.5)

    return total

end

function renormalon_MS_func(; z::Float64, Q::Float64, μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64, Ω1::Float64)

    total = Ω1/(2*Q*(z*(1-z))^1.5)

    return total

end

function EEC_MSR_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    pert = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
    renormalon = renormalon_func(z=z, Q=Q, μ_ren=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf, Ω1=0.0) 
    return pert+renormalon
    
end

function EEC_MSR_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = EEC_MSR_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end

function test_func(; μ_ren::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    R1 = μ_ren
    μ = μ_ren

    γR0 = γR0_func(nf)
    γR1 = γR1_func(nf)

    R_max = [R1]
    R_min = [R0]

    f1(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)/(4π)
    integral1, error1 = hcubature(f1, R_min, R_max)

    f2(x) = (alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)/(4π))^2
    integral2, error2 = hcubature(f2, R_min, R_max)

    order1 = integral1
    order2 = integral2

    if order == 1 
        total = order1
    end
    if order == 2 
        total = order2
    end

    return total

end
#a=Ω_Integral_func(μ_ren=91.2, μ_ini=91.2, αs_ini=0.118, order=2, nf=5)
#b=ΣR_func(μ_ren=91.2, μ_ini=91.2, αs_ini=0.118, order=2, nf=5)
#println(b-a)