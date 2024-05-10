# Non-resummed EEC in back to back limit https://arxiv.org/pdf/2012.07859

include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")

function L0_func(z)
    return 1/(1-z)
end

function L1_func(z)
    return log(1-z)/(1-z)
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

function SCET_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order, nf=nf)

    sigma1 = sigma1_func(z)

    if order == 1 
        total = (αs/(4π))*sigma1
    end

    return total

end

function SCET_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end