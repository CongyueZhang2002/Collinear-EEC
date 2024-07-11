# Non-resummed EEC in back to back limit https://arxiv.org/pdf/2012.07859

include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")

function L0_func(z)
    return 1/(1-z)
end

function L1_func(z)
    return log(1-z)/(1-z)
end

function L2_func(z)
    return (log(1-z))^2/(1-z)
end

function L3_func(z)
    return (log(1-z))^3/(1-z)
end

function L4_func(z)
    return (log(1-z))^4/(1-z)
end

function sigma1_func(z) 

    L0 = L0_func(z)
    L1 = L1_func(z)

    value = CF*(
        - 2*L1 
        - 3*L0
    )

    return value
end

function sigma2_func(z,nf)

    L0 = L0_func(z)
    L1 = L1_func(z)
    L2 = L2_func(z)
    L3 = L3_func(z)

    value = CF*(
        4 * CF * L3 
        + L2 * (18*CF + (22/3)*CA - (4/3)*nf) 
        + L1 * (CF*(34 + 8z2) + CA*(-35/9 + 4z2) + (2/9)*nf) 
        + L0 * (CF*(45/2 + 24*z2 - 8*z3) + CA*(-35/2 + 22*z2 + 12*z3) + nf*(3 - 4z2))
    ) 
    return value
end

function SCET_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    sigma1 = sigma1_func(z)
    sigma2 = sigma2_func(z,nf)

    if order == 1 
        total = (αs/(4π))*sigma1
    elseif order == 2 
        total = (αs/(4π))*sigma1 + (αs/(4π))^2*sigma2
    end

    return total

end

function SCET_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end