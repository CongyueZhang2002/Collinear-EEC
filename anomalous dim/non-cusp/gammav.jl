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

function γν1_func(Lb,nf)

    γν0 = γν1_func(nf)

    Γ0 = 4*Γ1_func(nf)

    order1 = -2*Lb*Γ0 + γν0

    return order1
end

function γν2_func(Lb,nf)

    γν0 = γν1_func(nf)
    γν1 = γν2_func(nf)

    β0 = 4*β0_func(nf)

    Γ0 = 4*Γ1_func(nf)
    Γ1 = 4^2*Γ2_func(nf)

    order2 = (
        - Lb^2*Γ0*β0 
        + Lb*(β0*γν0 - 2*Γ1) 
        + γν1
    )

    return order2
end

function γν3_func(Lb,nf)

    γν0 = γν1_func(nf)
    γν1 = γν2_func(nf)
    γν2 = γν3_func(nf)

    β0 = 4*β0_func(nf)
    β1 = 4^2*β1_func(nf)    

    Γ0 = 4*Γ1_func(nf)
    Γ1 = 4^2*Γ2_func(nf)
    Γ2 = 4^3*Γ3_func(nf)

    order3 = (
        - 2/3*Lb^3 * Γ0*β0^2 
        + Lb^2 * (β0^2*γν0 - 2*Γ1*β0 - Γ0*β1) 
        + Lb * (2β0*γν1 + β1*γν0 - 2*Γ2)
        + γν2
    )
    return order3
end

function γν_func(; b::Float64, μ0::Float64, αs::Float64, order::Int64, nf::Int64)

    Lb = log((b*μ0/b0)^2)

    γν1 = γν1_func(Lb, nf)
    γν2 = γν2_func(Lb, nf)
    γν3 = γν3_func(Lb, nf)  

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

function γν_final(; b::Float64, μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64, bmax::Float64)

    bstar = b/(1+(b/bmax)^2)^0.5

    μ0 = b0/bstar

    αs_μ0 = alpha_s_func(μf=μ0, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

    γν_FO = γν_func(b=b, μ0=μ0, αs=αs_μ0, order=order, nf=nf)

    μ_max = [μ]
    μ_min = [μ0]

    f(x) = -0.25/x[1]*Γ_final(; μ=x[1], μ_ini=μ_ini, αs_ini=αs_ini, order=order+1, nf=nf)

    integral, error = hcubature(f, μ_min, μ_max)

    total = integral + γν_FO

    return total

end