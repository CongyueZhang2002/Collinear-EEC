using Distributed
using HCubature
using Plots

include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\..\\strong coupling\\constants.jl")

include("gammaB.jl")
include("gammaS.jl")
include("..\\cusp\\cusp.jl")

@everywhere function γH_final(; Q::Float64, μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
    Γ = Γ_func(αs=αs, order=order+1, nf=nf)

    γS = -γS_func(αs=αs, order=order, nf=nf)
    γJ = γB_func(αs=αs, order=order, nf=nf) + γS_func(αs=αs, order=order, nf=nf)

    γH_f = 4*Γ*log(Q/μ) - 2γJ - γS
    
    return γH_f
end

@everywhere function log_αs(; αs_μi::Float64, αs_μf::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)

    β_1(x) = 1/(-2*x[1]*(x[1]/(4π)*β0))
    β_2(x) = 1/(-2*x[1]*(x[1]/(4π)*β0 + (x[1]/(4π))^2*β1))
    β_3(x) = 1/(-2*x[1]*(x[1]/(4π)*β0 + (x[1]/(4π))^2*β1 + (x[1]/(4π))^3*β2))
    β_4(x) = 1/(-2*x[1]*(x[1]/(4π)*β0 + (x[1]/(4π))^2*β1 + (x[1]/(4π))^3*β2 + (x[1]/(4π))^4*β3))

    if order == 1 
        integral, error = hcubature(β_1, [αs_μi], [αs_μf])
    end
    if order == 2 
        integral, error = hcubature(β_2, [αs_μi], [αs_μf])
    end
    if order == 3 
        integral, error = hcubature(β_3, [αs_μi], [αs_μf])
    end
    if order == 4 
        integral, error = hcubature(β_4, [αs_μi], [αs_μf])
    end
      
    return integral
end

@everywhere function γH_analytic(; Q::Float64, αs_μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs_Q = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

    Γ = Γ_func(αs=αs_μ, order=order+1, nf=nf)
    γS = -γS_func(αs=αs_μ, order=order, nf=nf)
    γJ = γB_func(αs=αs_μ, order=order, nf=nf) + γS_func(αs=αs_μ, order=order, nf=nf)

    log_term = log_αs(αs_μi=αs_μ, αs_μf=αs_Q, order=order, nf=nf)

    γH_f = 4*Γ*log_term - 2γJ - γS
    
    return γH_f
end


