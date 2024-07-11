# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

using Distributed
using HCubature

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")
include("..\\non-cusp\\Integration.jl")

# rapidity ν

@everywhere function γν0_func(nf)

    γν0 = 0

    return γν0
end

@everywhere function γν1_func(nf) 

    #β0 = β0_func(nf) # different convention 1/4π vs 1/π
    #γν1 = 2CF * (CA*(-64/9 + 28z3) - β0*56/9)

    γν1 = 4*(7.46333 + 2.76543*nf)

    return γν1
end

@everywhere function γν2_func(nf)

    #β0 = β0_func(nf)
    #β1 = β1_func(nf)
    #γν2 = 2CF * (
    #    + CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
    #    + CA*β0 * (3865/54 + 412/27*z2 + 220/9*z3 - 50z4) 
    #    + β0^2 * (-464/81 - 8z3) 
    #    + β1 * (-1711/54 + 152/9*z3 + 8z4)
    #)

    γν2 = 4*(70.068 + 77.1286*nf - 4.54662*nf^2)

    return γν2
end

@everywhere function γν3_func(nf)

    γν3 = 4*(- 350.8 + 2428*nf - 378.3*nf^2 + 8.072*nf^3)

    return γν3
end

@everywhere function γν_raw(; αs::Float64, order::Int64, nf::Int64)

    γν0 = γν0_func(nf)
    γν1 = γν1_func(nf)
    γν2 = γν2_func(nf)
    γν3 = γν3_func(nf)

    order1 = αs/(4π) * γν0
    order2 = (αs/(4π))^2 * γν1
    order3 = (αs/(4π))^3 * γν2
    order4 = (αs/(4π))^4 * γν3

    if order == 1 
        total = order1
    elseif order == 2  
        total = order1 + order2
    elseif order == 3 
        total = order1 + order2 + order3
    elseif order == 4 
        total = order1 + order2 + order3 + order4 
    end

    return total
end

@everywhere function γν0_final(Lb,nf)

    γν0 = γν0_func(nf)

    Γ0 = Γ0_func(nf)

    order1 = -2*Lb*Γ0 + γν0

    return order1
end

@everywhere function γν1_final(Lb,nf)

    γν0 = γν0_func(nf)
    γν1 = γν1_func(nf)

    β0 = β0_func(nf)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)

    order2 = (
        - Lb^2*Γ0*β0 
        + Lb*(β0*γν0 - 2*Γ1) 
        + γν1
    )

    return order2
end

@everywhere function γν2_final(Lb,nf)

    γν0 = γν0_func(nf)
    γν1 = γν1_func(nf)
    γν2 = γν2_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)    

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)

    order3 = (
        - 2/3*Lb^3 * Γ0*β0^2 
        + Lb^2 * (β0^2*γν0 - 2*Γ1*β0 - Γ0*β1) 
        + Lb * (2β0*γν1 + β1*γν0 - 2*Γ2)
        + γν2
    )
    return order3
end

@everywhere function γν3_final(Lb,nf)

    γν0 = γν0_func(nf)
    γν1 = γν1_func(nf)
    γν2 = γν2_func(nf)
    γν3 = γν3_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)    
    β2 = β2_func(nf)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)

    order4 = (
        - 1/2*Lb^4 * Γ0*β0^3 
        + Lb^3 * (β0^3*γν0 - 2*β0^2*Γ1 -5/3*β1*β0*Γ0)
        + Lb^2 * (3*β1^2*γν1 + 5/2*β1*β0*γν0 - β2*Γ0 - 2*Γ1*β1 - 3*β0*Γ2) 
        + Lb * (β2*γν0 + 2*β1*γν1 + 3*β0*γν2 - 2*Γ3)
        + γν3
    )
    return order4
end

@everywhere function γν_func(; b::Float64, μ0::Float64, αs::Float64, order::Int64, nf::Int64, bmax::Float64)

    #bstar = b/(1+(b/bmax)^2)^0.5

    Lb = log((b*μ0/b0)^2)

    γν0 = γν0_final(Lb, nf)
    γν1 = γν1_final(Lb, nf)
    γν2 = γν2_final(Lb, nf)
    γν3 = γν3_final(Lb, nf)

    order1 = αs/(4π) * γν0
    order2 = (αs/(4π))^2 * γν1
    order3 = (αs/(4π))^3 * γν2
    order4 = (αs/(4π))^4 * γν3

    if order == 1 
        total = order1
    elseif order == 2  
        total = order1 + order2
    elseif order == 3 
        total = order1 + order2 + order3
    elseif order == 4 
        total = order1 + order2 + order3 + order4 
    end

    return total

end

#numerical integration

@everywhere function γν_final(; b::Float64, μB::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64, bmax::Float64)

    bstar = b/(1+(b/bmax)^2)^0.5

    μ0 = b0/bstar

    αs_μ0 = alpha_s_func(μf=μ0, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    γν_FO = γν_func(b=b, μ0=μ0, αs=αs_μ0, order=order, nf=nf, bmax=bmax)

    μ_max = [μB]
    μ_min = [μ0] 

    f(x) = -4/x[1]*Γ_final(; μ=x[1], μ_ini=μ_ini, αs_ini=αs_ini, order=order+1, nf=nf)

    integral, error = hcubature(f, μ_min, μ_max)

    total = integral + γν_FO

    return total

end