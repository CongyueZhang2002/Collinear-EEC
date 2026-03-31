# Non-resummed EEC in back to back limit https://arxiv.org/pdf/2012.07859

#include("..\\strong coupling\\constants.jl")
#include("..\\strong coupling\\alpha_ss.jl")

function FO_func(; z::Float64, Q::Float64)

    order = nloops_FO

    αs = alpha_s_func(μren_ratio*Q)

    Lz = log(z)
    LQ = log(μren_ratio^2)

    A = (2)/2
    B = (-11.5333333333333333333333333333*Lz + 81.4809061031774183115800436626)/4
    C = (45.1488888888888888888888888889*Lz^2 -1037.73139012054533703493186893*Lz + 2871.36018216039823956603951601)/8

    β0 = β0_func(nf)
    β1 = β1_func(nf)

    order1 = (αs/(2π))*A
    order2 = (αs/(2π))^2*(
          B 
        + 1/2*β0*LQ*A
        )
    order3 = (αs/(2π))^3*(
          C 
        + β0*LQ*B
        + (1/4*β1*LQ + 1/4*β0^2*LQ^2)*A
        )  

    if order == 0 
        total = order1
    elseif order == 1 
        total = order1 + order2
    elseif order == 2 
        total = order1 + order2 + order3
    end

    #if order == 0
    #    value = (2*as)/z
    #elseif order == 1 
    #    value = (2*as + as^2*(-11.5333333333333333333333333333*Lz + 81.4809061031774183115800436626))/z
    #elseif order == 2
    #    value = (
    #      2*as 
    #    + as^2*(-11.5333333333333333333333333333*Lz + 81.4809061031774183115800436626) 
    #    + as^3*(45.1488888888888888888888888889*Lz^2 -1037.73139012054533703493186893*Lz + 2871.36018216039823956603951601)
    #    )/z
    #end

    return total/z
end

function FO_sigma_χ(; χ::Float64, Q::Float64)

    z = 0.5*(1-cos(χ))

    part = FO_func(z=z, Q=Q)

    total = 0.5*sin(χ)*part

    return total

end