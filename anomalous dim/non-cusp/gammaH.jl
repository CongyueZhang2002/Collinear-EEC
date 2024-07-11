using Distributed
using HCubature
using Plots

include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\..\\strong coupling\\constants.jl")

include("gammaB.jl")
include("gammaS.jl")
include("..\\cusp\\cusp.jl")

@everywhere function γH0_tilde(nf)

    γS0_t = γS0_tilde(nf)
    γJ0_t = γB0_tilde(nf) 

    γH0_t = - 2γJ0_t - γS0_t
    
    return γH0_t
end

@everywhere function γH1_tilde(nf)

    γS1_t = γS1_tilde(nf)
    γJ1_t = γB1_tilde(nf) 

    γH1_t = - 2γJ1_t - γS1_t
    
    return γH1_t
end

@everywhere function γH2_tilde(nf)

    γS2_t = γS2_tilde(nf)
    γJ2_t = γB2_tilde(nf) 

    γH2_t = - 2γJ2_t - γS2_t
    
    return γH2_t
end

@everywhere function γH3_tilde(nf)

    γS3_t = γS3_tilde(nf)
    γJ3_t = γB3_tilde(nf) 

    γH3_t = - 2γJ3_t - γS3_t
    
    return γH3_t
end

@everywhere function γH_tilde(; αs::Float64, order::Int64, nf::Int64)

    γS_t = γS_tilde(αs=αs, order=order, nf=nf)
    γJ_t = γB_tilde(αs=αs, order=order, nf=nf) 

    γH_t = - 2γJ_t - γS_t
    
    return γH_t
end

#numerical integration

#@everywhere function γH_final(; Q::Float64, μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
#
#    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
#    Γ = Γ_func(αs=αs, order=order+1, nf=nf)
#
#    γS = -γS_func(αs=αs, order=order, nf=nf)
#    γJ = γB_func(αs=αs, order=order, nf=nf) + γS_func(αs=αs, order=order, nf=nf)
#
#    γH_f = 4*Γ*log(Q/μ) - 2γJ - γS
#    
#    return γH_f
#end


