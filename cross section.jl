using Distributed

@everywhere begin

using QuadGK
using HCubature
using SpecialFunctions
include("fixed order/analytic.jl")

function Integrand_func(; b::Float64, z::Float64, Q::Float64, Coeffs::Vector{Float64})

    k = sqrt(z)*Q
    bstar = b/sqrt(1+b^2/b0^2)

    Lb = log((μren_ratio*Q/(μJ_ratio*b0/bstar))^2)

    Lb_vec = [1.0, Lb, Lb^2, Lb^3, Lb^4, Lb^5, Lb^6, Lb^7, Lb^8, Lb^9]
    Σb = sum(Lb_vec .* Coeffs)

    NP_Sudakov = NP(b)
    value = Q^2/2*b*besselj0(b*k)*(Σb*NP_Sudakov)

    return value
end

function sigma_z(; z::Float64, Q::Float64)

    as = alpha_s_func(μren_ratio*Q)/(4π)
    
    Coeffs = Σb_func(as)

    integrand(b) = Integrand_func(b=b, z=z, Q=Q, Coeffs=Coeffs)

    total, error = quadgk(integrand, 0.0001, 5.0, rtol=rtol)

    return total
end

function R_ratio(Q)
    nf = 5
    αs = alpha_s_func(Q)
    r = 1 + αs/π + (αs/(4π))^2*(CF*CA*(123/2-44*z3) + CF*TF*nf*(-22+16*z3) - CF^2*3/2)
    return r
end

function sigma_χ(; χ::Float64, Q::Float64)

    z = 0.5*(1-cos(χ))

    part = sigma_z(z=z, Q=Q) 

    total = 0.5*sin(χ)*part

    r = R_ratio(Q)

    return total/r
end

end