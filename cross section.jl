using Distributed

@everywhere begin

using QuadGK
using HCubature
using SpecialFunctions
include("fixed order/analytic.jl")
#include("strong coupling/constants.jl") 
#include("strong coupling/alpha_s.jl")
#include("NP.jl")

function Integrand_func(; b::Float64, z::Float64, Q::Float64, 
    as::Float64, order::Int64, nf::Int64, 
    params::Vector{Float64}, μJ_ratio::Float64, μH_ratio::Float64,
    Coeffs::Vector{Float64})

    μ = Q

    k = sqrt(z)*Q
    bstar = b/sqrt(1+b^2/b0^2)

    ζi = (b0/bstar)^2
    ζf = μ^2

    Lb = log(bstar^2*μ^2/(μJ_ratio^2*b0^2))

    Lb_vec = [1, Lb, Lb^2, Lb^3, Lb^4, Lb^5, Lb^6, Lb^7, Lb^8, Lb^9]
    Σb = sum(Lb_vec .* Coeffs)

    NP_Sudakov = NP(z=z, b=b, params=params, ζi=ζi, ζf=ζf)
    value = Q^2/2*b*besselj0(b*k)*(Σb*NP_Sudakov)

    return value
end

function sigma_z(; z::Float64, Q::Float64, 
    μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64, 
    params::Vector{Float64}=[0.0], μJ_ratio::Float64, μH_ratio::Float64)

    as = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)/(4π)
    
    Coeffs = Σb_func(as=as, order=order, nf=nf, μJ_ratio=μJ_ratio, μH_ratio=μH_ratio)

    integrand(b) = Integrand_func(b=b, z=z, Q=Q, 
        as=as, order=order, nf=nf, 
        params=params, μJ_ratio=μJ_ratio, μH_ratio=μH_ratio,
        Coeffs=Coeffs
    )

    total, error = quadgk(integrand, 0.0001, 5.0, rtol=1e-5)

    return total
end

function R_ratio(; Q::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64)
    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=3, nf=nf)
    r = 1 + αs/π + (αs/(4π))^2*(CF*CA*(123/2-44*z3) + CF*TF*nf*(-22+16*z3) - CF^2*3/2)
    return r
end

function sigma_χ(; χ::Float64, Q::Float64, 
    μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64, 
    params::Vector{Float64}=[0.0], μJ_ratio::Float64, μH_ratio::Float64)

    z = 0.5*(1-cos(χ))

    part = sigma_z(z=z, Q=Q, 
    μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf, params=params
    , μJ_ratio=μJ_ratio, μH_ratio=μH_ratio) 

    total = 0.5*sin(χ)*part

    r = R_ratio(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf)

    return total/r
end

end