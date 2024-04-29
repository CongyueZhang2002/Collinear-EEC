include("gammaB.jl")
include("gammaS.jl")
include("..\\cusp\\cusp.jl")

function γH_f_func(; Q::Float64, μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
    Γ = Γ_func(αs=αs, order=order+1, nf=nf)

    γS = -γS_func(αs=αs, order=order, nf=nf)
    γJ = γB_func(αs=αs, order=order, nf=nf) + γS_func(αs=αs, order=order, nf=nf)

    γH_f = 4*Γ*log(Q/μ) - 2γJ - γS
    
    return γH_f
end