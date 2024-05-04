# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

using Distributed
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")

# Soft

@everywhere function γS0_func(nf)

    γS0 = 0

    return γS0
end

@everywhere function γS1_func(nf) 

    β0 = β0_func(nf) 

    γS1 = 2CF * (CA*(-64/9 + 28z3) + β0*(-56/9 + 2z2))

    return γS1
end

@everywhere function γS2_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)

    γS2 = 2CF * (
          CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
        + CA*β0 * (4697/54 + 484/27*z2 + 220/9*z3 - 112z4) 
        + β0^2 * (520/81 + 10/3*z2 - 28/3*z3) 
        + β1 * (-1711/54 + 2z2 + 152/9*z3 + 8z4)
    )

    return γS2
end

@everywhere function γS_func(; αs::Float64, order::Int64, nf::Int64)

    γS0 = γS0_func(nf) 
    γS1 = γS1_func(nf)
    γS2 = γS2_func(nf) 

    order1 = αs/(4π) * γS0
    order2 = (αs/(4π))^2 * γS1
    order3 = (αs/(4π))^3 * γS2

    if order == 1 
        total = order1
    end
    if order == 2  
        total = order1 + order2
    end
    if order == 3 
        total = order1 + order2 + order3
    end

    return total

end

# final γS in integration
@everywhere function γS_final(; μ::Float64, ν::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
    Γ = Γ_func(αs=αs, order=order+1, nf=nf)
    γS = -γS_func(αs=αs, order=order, nf=nf)

    γS_f = 4*Γ*log(μ/ν) + γS
    
    return γS_f
end