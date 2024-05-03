# cusp anomalous dimension up to 4 loops from https://arxiv.org/pdf/1911.10174.pdf

using Distributed
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")

@everywhere function Γ1_func(nf)
    Γ1 = CF
    return Γ1
end

@everywhere function Γ2_func(nf)
    Γ2 = CF * (1.03864*CA - 0.555556*nf*TF)
    return Γ2
end

@everywhere function Γ3_func(nf)
    Γ3 = CF * (1.52982*CA^2 - 1.45614*CA*nf*TF + 0.0562236*CF*nf*TF - 0.0370370*nf^2*TF^2)
    return Γ3
end

@everywhere function Γ4_func(nf)
    Γ4 = CF * (
        2.38379*CA^3 - 3.44271*CA^2*nf*TF + 0.303089*CA*CF*nf*TF - 0.242621*CF^2*nf*TF 
        + 0.911990*CA*nf^2*TF^2 - 0.333037*CF*nf^2*TF^2 + 0.0766956*nf^3*TF^3
    ) - dFdAdNA*NA/NF*1.97915 - nf*dFdFdNA*NA/NF*0.483964
    return Γ4
end

@everywhere function Γ_func(; αs::Float64, order::Int64, nf::Int64)

    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)    

    order1 = αs/π * Γ1
    order2 = (αs/π)^2 * Γ2
    order3 = (αs/π)^3 * Γ3
    order4 = (αs/π)^4 * Γ4

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

    return total

end

@everywhere function Γ_final(; μ::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order, nf=nf)

    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)    

    order1 = αs/π * Γ1
    order2 = (αs/π)^2 * Γ2
    order3 = (αs/π)^3 * Γ3
    order4 = (αs/π)^4 * Γ4

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

    return total

end