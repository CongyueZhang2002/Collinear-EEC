using Distributed

@everywhere begin

using QuadGK
using SpecialFunctions
include("fixed order/analytic.jl")

# Interpolations

function Integrand_func(; b::Float64, z::Float64, Q::Float64, as::Float64)

    k = sqrt(z)*Q
    bstar = b/sqrt(1+b^2/b0^2)

    Lb = log((μren_ratio*Q/(μJ_ratio*b0/bstar))^2)

    a = params

    init_splitting_globals(a1b=a[2]*b, a2=a[3])
    init_hard_globals(a1b=a[2]*b, a2=a[3])

    Σb_quark = Σb_quark_func(Lb=Lb, a_s=as, order=2, nf=5)   
    Σb_gluon = Σb_gluon_func(Lb=Lb, a_s=as, order=2, nf=5)

    value = Q^2/2*b*besselj0(b*k)*(Σb_quark + Σb_gluon)

    return value
end

function sigma_z(; z::Float64, Q::Float64)

    as = alpha_s_func(μren_ratio*Q)/(4π)

    integrand(b) = Integrand_func(b=b, z=z, Q=Q, as=as)

    #println("One Integration")

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