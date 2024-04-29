# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

using HCubature

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")

# rapidity ν

function γν1_func(nf)

    γν1 = 0

    return γν1
end

function γν2_func(nf) 

    β0 = 4*β0_func(nf) # different convention 1/4π vs 1/π

    γν2 = 2CF * (CA*(-64/9 + 28z3) - β0*56/9)

    return γν2
end

function γν3_func(nf)

    β0 = 4*β0_func(nf)
    β1 = 4^2*β1_func(nf)

    γν3 = 2CF * (
        + CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
        + CA*β0 * (3865/54 + 412/27*z2 + 220/9*z3 - 50z4) 
        + β0^2 * (-464/81 - 8z3) 
        + β1 * (-1711/54 + 152/9*z3 + 8z4)
    )

    return γν3
end


function γν_func(; αs::Float64, order::Int64, nf::Int64)

    γν1 = γν1_func(nf)
    γν2 = γν2_func(nf)
    γν3 = γν3_func(nf)  

    order1 = αs/(4π) * γν1
    order2 = (αs/(4π))^2 * γν2
    order3 = (αs/(4π))^3 * γν3

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

function Γv_func(; μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

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

function γν_f_func(; b::Float64, μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    μ0 = b0/b

    integrand(μ') = 1/μ'*Γv_func(; μ=μ', αs_ini=αs_ini, μ_ini=μ_ini, order=order+1, nf=nf)
    xmax = [μ]
    xmin = [μ0]
    integral, error = hcubature(integrand, xmin, xmax)

    αs = alpha_s_func(μf=μ0, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
    γν = γν_func(αs=αs, order=order, nf=nf)

    γν_f = -4*integral + γν
    
    return γν_f
end