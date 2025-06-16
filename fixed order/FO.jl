# Non-resummed EEC in back to back limit https://arxiv.org/pdf/2012.07859

#include("..\\strong coupling\\constants.jl")
#include("..\\strong coupling\\alpha_ss.jl")

function FO_func(; z::Float64, Q::Float64)

    order = nloops_FO

    as = alpha_s_func(μren_ratio*Q)/(4π)

    Lz = log(z/μren_ratio^2)

    if order == 0
        value = (2*as)/z
    elseif order == 1 
        value = (2*as + as^2*(-11.5333333333333333333333333333*Lz + 81.4809061031774183115800436626))/z
    elseif order == 2
        value = (
          2*as 
        + as^2*(-11.5333333333333333333333333333*Lz + 81.4809061031774183115800436626) 
        + as^3*(45.1488888888888888888888888889*Lz^2 -1037.73139012054533703493186893*Lz + 2871.36018216039823956603951601)
        )/z
    end

    return value
end

function FO_sigma_χ(; χ::Float64, Q::Float64)

    z = 0.5*(1-cos(χ))

    part = FO_func(z=z, Q=Q)

    total = 0.5*sin(χ)*part

    return total

end