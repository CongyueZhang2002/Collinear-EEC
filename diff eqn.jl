using DifferentialEquations
using ForwardDiff
include("strong coupling\\constants.jl")
include("strong coupling\\alpha_s.jl") 

function coupled_system_old!(dv, v, params, μ)
    jq, jg = v  # unpack state variables

    nf = 5
    as = alpha_s_resum(μf=μ, μi=91.2, αs=0.118, order=2, nf=nf)/(4*pi)

    γ_qq = 25/6*CF
    γ_gq = -7/6*CF 
    γ_qg = -7/15*nf
    γ_gg = 14/5*CA+2/3*nf

    djq = -2/μ*as*(γ_qq*jq+γ_gq*jg)
    djg = -2/μ*as*(γ_qg*jq+γ_gg*jg)

    dv[1] = djq
    dv[2] = djg
end

function collinear_jet_old(; b::Float64, Q::Float64)

    dv0 = [1.0,1.0]

    μspan = (b0/b, Q)

    problem = ODEProblem(coupled_system_old!, dv0, μspan)
    solution = solve(problem)

    return solution.u[end][1]
end

function Integrand_old(; b::Float64, z::Float64, Q::Float64, 
    params::Vector{Float64})

    k = sqrt(z)*Q
    bstar = b/sqrt(1+b^2/b0^2)

    ζi = (b0/bstar)^2
    ζf = Q^2

    Σb = 1/2*collinear_jet_old(b=bstar, Q=Q)

    NP_Sudakov = NP(b=b, params=params, ζi=ζi, ζf=ζf)
    value = Q^2/2*b*besselj0(b*k)*Σb*NP_Sudakov

    return value
end

function dσdz_old(; z::Float64, Q::Float64, params::Vector{Float64}=[0.0]) 
             
    integrand(b) = Integrand_old(b=b[1], z=z, Q=Q, params=params)

    total, error = hcubature(integrand, [0.0001], [10.0])

    return total
end


function dσdz_plot(; z::Float64, Q::Float64) 

    value = dσdz_old(z=z, Q=Q) 
    return value
end

#println(dσdz_old(z=0.1, Q=91.2))