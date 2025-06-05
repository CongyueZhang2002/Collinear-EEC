# Non-resummed EEC in back to back limit https://arxiv.org/pdf/2012.07859

#include("..\\strong coupling\\constants.jl")
#include("..\\strong coupling\\alpha_ss.jl")

function SCET_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    as = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)/(4π)

    if order == 0
        value = (2*as)/z
    elseif order == 1 
        value = (2*as + as^2*(-11.5333333333333333333333333333*log(z) + 81.4809061031774183115800436626))/z
    elseif order == 2
        value = (
        2*as 
        + as^2*(-11.5333333333333333333333333333*log(z) + 81.4809061031774183115800436626) 
        + as^3*(45.1488888888888888888888888889*log(z)^2 -1037.73139012054533703493186893*log(z) +2871.36018216039823956603951601)
        )/z
    end

    return value
end

function SCET_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end