# cusp anomalous dimension up to 4 loops from https://arxiv.org/pdf/1911.10174.pdf
# expand in order of (αs/4π)^n
# https://arxiv.org/pdf/2205.02242 5 loop apporximation of cusp
using Distributed
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")

@everywhere function Γ0_func(nf)
    Γ0 = CF
    return 4*Γ0
end

@everywhere function Γ1_func(nf)
    Γ1 = CF * (1.03864*CA - 0.555556*nf*TF)
    return 4^2*Γ1
end

@everywhere function Γ2_func(nf)
    Γ2 = CF * (1.52982*CA^2 - 1.45614*CA*nf*TF + 0.0562236*CF*nf*TF - 0.0370370*nf^2*TF^2)
    return 4^3*Γ2
end

@everywhere function Γ3_func(nf)
    Γ3 = CF * (
        2.38379*CA^3 - 3.44271*CA^2*nf*TF + 0.303089*CA*CF*nf*TF - 0.242621*CF^2*nf*TF 
        + 0.911990*CA*nf^2*TF^2 - 0.333037*CF*nf^2*TF^2 + 0.0766956*nf^3*TF^3
    ) - dFdAdNA*NA/NF*1.97915 - nf*dFdFdNA*NA/NF*0.483964
    return 4^4*Γ3
end

@everywhere function Γ4_func(nf)
    return 0.21 # 0.21+-0.17 @ nf=5
end

@everywhere function Γ_func(; αs::Float64, order::Int64, nf::Int64)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)    

    order1 = αs/(4π) * Γ0
    order2 = (αs/(4π))^2 * Γ1
    order3 = (αs/(4π))^3 * Γ2
    order4 = (αs/(4π))^4 * Γ3
    order5 = (αs/(4π))^5 * Γ4

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

    return total

end

@everywhere function Γ_final(; μ::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)    

    order1 = αs/(4π) * Γ0
    order2 = (αs/(4π))^2 * Γ1
    order3 = (αs/(4π))^3 * Γ2
    order4 = (αs/(4π))^4 * Γ3
    order5 = (αs/(4π))^5 * Γ4

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

    return total

end