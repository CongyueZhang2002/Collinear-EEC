# https://arxiv.org/pdf/2403.04077 equation 10

include("..\\strong coupling\\alpha_s.jl")

#function A1_func(z)
#    value = 4/3 * (
#        - 1/2*log(1-z)/(1-z) 
#        - 3/4*(1-z) 
#        + 3/8*(1/z)
#        + 1/(8z^5) * (
#            4*(- z^4 - z^3 + 3z^2 - 15z + 9)*log(1-z)
#            - 9z^4 - 6z^3 - 42z^2 + 36z
#        )
#    )
#    return value
#end

function A1_func(z)
    value = CF * (3-2z)/(4*(1-z)*z^5)*(
        3z*(2-3z)+2*(2z^2-6z+3)*log(1-z)
    )
    return value
end

function perturbation_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order, nf=nf)

    A1 = A1_func(z)

    if order == 1 
        total = (αs/π)*A1
    end

    return total

end

function perturbation_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end
