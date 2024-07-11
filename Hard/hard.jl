# hard function up to 3 loops from https://arxiv.org/pdf/1004.3653

using Distributed
include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")

function γ0H_func(nf)

    γ0H = -6*CF

    return γ0H
end

function γ1H_func(nf)

    γ1H=(
          CF^2 * (-3 + 4*π^2 - 48*z3) 
        + CF*CA * (-961/27 - 11*π^2/3 + 52*z3) 
        + CF*TF*nf * (260/27 + 4*π^2/3)
    )

    return γ1H
end

function γ2H_func(nf)

    γ2H= (-1.856*nf^2 + 259.3*nf - 1499)

    return γ2H
end

function c1H_func(nf) 

    c1H = CF*(-16 + 7*π^2/3)

    return c1H
end

function c2H_func(nf) 

    c2H = (
          CF^2 * (511/4 - 83*π^2/3 + 67*π^4/30 - 60*z3) 
        + CF*CA * (-51157/324 + 1061*π^2/54 - 8*π^4/45 + 626*z3/9) 
        + CF*TF*nf * (4085/81 - 182*π^2/27 + 8*z3/9)
    )

    return c2H    
end

@everywhere function H1_func(L,nf) 

    γ0H = γ0H_func(nf)

    c1H = c1H_func(nf)

    Γ0 = Γ0_func(nf)

    H1 = - 1/2*Γ0*L^2 - γ0H*L + c1H

    return H1
end

@everywhere function H2_func(L,nf)

    γ0H = γ0H_func(nf)
    γ1H = γ1H_func(nf)

    c1H = c1H_func(nf)
    c2H = c2H_func(nf)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)

    β0 = β0_func(nf)

    H2 = (
          1/8*Γ0^2 * L^4 
        + (β0*Γ0/6 + γ0H*Γ0/2) * L^3 
        + (γ0H^2/2 + β0*γ0H/2 - Γ1/2) * L^2 
        - γ1H * L 
        + c1H * (-L^2*Γ0/2 + L*(-β0 - γ0H))
        + c2H
    )

    return H2
end

@everywhere function H3_func(L,nf) 

    γ0H = γ0H_func(nf)
    γ1H = γ1H_func(nf)
    γ2H = γ2H_func(nf)

    c1H = c1H_func(nf)
    c2H = c2H_func(nf)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)

    H3 = (
        - 1/48*Γ0^3 * L^6 
        + (-1/12*β0*Γ0^2 - 1/8*γ0H*Γ0^2) * L^5 
        + (-1/12*Γ0*β0^2 - 5/12*γ0H*Γ0*β0 - 1/4*γ0H^2*Γ0 + Γ0*Γ1/4) * L^4 
        + (-γ0H^3/6 - 1/2*β0*γ0H^2 - 1/3*β0^2*γ0H + Γ1*γ0H/2 + β1*Γ0/6 + γ1H*Γ0/2 + β0*Γ1/3) * L^3 
        + (β1*γ0H/2 + γ1H*γ0H + β0*γ1H - Γ2/2) * L^2 
        - γ2H * L 
        + c1H * (
              1/8*Γ0^2*L^4 
            + (2*β0*Γ0/3 + γ0H*Γ0/2)*L^3 
            + (β0^2 + 3*γ0H*β0/2 + γ0H^2/2 - Γ1/2) * L^2 
            + (-β1 - γ1H) * L
        ) 
        + c2H * (L*(-2*β0 - γ0H) - L^2*Γ0/2) 
        # + c3H not yet known but doesn't really matter
    )

    return H3
end

@everywhere function H_func(; μH::Float64, Q::Float64, αs::Float64, order::Int64, nf::Int64)

    L = log(Q^2/μH^2)

    H1 = H1_func(L,nf)
    H2 = H2_func(L,nf)
    H3 = H3_func(L,nf)

    if order == 1 
        total = 1 + (αs/4π)*H1
    elseif order == 2  
        total = 1 + (αs/4π)*H1 + (αs/4π)^2*H2
    elseif order == 3 
        total = 1 + (αs/4π)*H1 + (αs/4π)^2*H2 + (αs/4π)^3*H3
    end

    return total

end