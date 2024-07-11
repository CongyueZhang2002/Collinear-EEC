# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811
# 4 loop https://arxiv.org/pdf/2205.02249
# numerical https://arxiv.org/pdf/1912.12920 eqn(3.12); checked and proved γB, γS are correct from 1-4 loops

using Distributed
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")
include("gammaS.jl")

# Beam

@everywhere function γB0_func(nf)

    γB0 = 6CF

    return γB0
end

@everywhere function γB1_func(nf) 

    β0 = β0_func(nf) # different convention 1/4π vs 1/π

    γB1 = 2CF * (
          CA * (73/9 - 40z3) 
        + CF * (3/2 - 12z2 + 24z3) 
        + β0 * (121/18 + 2z2)
    )

    return γB1
end

@everywhere function γB2_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)

    γB2 = 2CF * (
          CA^2 * (52019/162 - 1682/27*z2 - 2056/9*z3 - 820/3*z4 + 176/3*z2*z3 + 232z5) 
        + CA*CF * (151/4 - 410/3*z2 + 844/3*z3 - 494/3*z4 + 16*z2*z3 + 120z5) 
        + CF^2 * (29/2 + 18z2 + 68z3 + 144z4 - 32z2*z3 - 240z5) 
        + CA*β0 * (-7739/54 + 650/27*z2 - 1276/9*z3 + 617/3*z4) 
        + β0^2 * (-3457/324 + 10/3*z2 + 16/3*z3) 
        + β1 * (1166/27 - 16/3*z2 + 52/9*z3 - 82/3*z4)
    )

    return γB2
end

@everywhere function γB0_tilde(nf)

    γB0 = γB0_func(nf)
    γS0 = γS0_func(nf)

    γB0_t = γB0 + γS0

    return γB0_t
end

@everywhere function γB1_tilde(nf)

    γB1 = γB1_func(nf)
    γS1 = γS1_func(nf)

    γB1_t = γB1 + γS1

    return γB1_t
end

@everywhere function γB2_tilde(nf)

    γB2 = γB2_func(nf)
    γS2 = γS2_func(nf)

    γB2_t = γB2 + γS2

    return γB2_t
end

@everywhere function γB3_tilde(nf)

    γB3_t = 2 * (
          bdFA*dFAdn
        + CA^3*CF * (-bdFA/24 - 371201/648 + 528*z3^2 + (8*z4 - 153670/81)*z3 - 11194/27*z4 + 6046/9*z6 + 11372/9*z5 + 472/3*z3*z2 + 504*z5*z2 + 4582/3*z2 - 2870*z7) 
        + nf*CA^2*CF * (-1/2*bNCACF2 - 1/4*bNCF3 - bdFF/48 - 16/3*z3^2 - 248/3*z5 - 137/9*z3 + 16186/27*z4 + (-584/9*z3 - 85175/162)*z2 - 144*z6 + 353/3) 
        + nf*CA*CF^2 * bNCACF2
        + nf*CF^3 * bNCF3
        + nf^2*CA*CF * (-320/9*z3 - 88/9*z5 - 80/9*z4 + (80/3*z3 + 3170/81)*z2 - 193/54)
        + nf^2*CF^2 * (-2104/27*z4 + 56/27*z3 + 368/9*z5 - 160/9*z3*z2 + 1244/27*z2 - 188/27) 
        + CA*CF^3 * (-2085/4 + 3220*z3^2 + (128*z4 - 3260)*z3 + 79297/18*z6 + 2167*z4 - 976*z5 + z2*(-1988/3*z3 + 2064*z5 + 1167) - 10920*z7) 
        + CF^2*CA^2 * (29639/36 - 7102/3*z3^2 + (129662/27 - 32*z4)*z3 + 5354/9*z5 - 60850/27*z4 - 5497/2*z6 + z2*(2096/9*z3 - 2104*z5 - 46771/27) + 8610*z7) 
        + bdFF*nf*dFFdn
        + CF*nf^3 * (-32/27*z4 + 32/81*z2 + 304/81*z3 - 131/81)
        + CF^4 * (-1152*z3^2 + 64*z4*z3 + 2004*z3 - 342*z4 + z2*(-120*z3 - 384*z5 - 450) - 2520*z5 - 2111*z6 + 5880*z7 + 4873/24)
    )

    return γB3_t
end

@everywhere function γB_tilde(; αs::Float64, order::Int64, nf::Int64)

    γB0_t = γB0_tilde(nf)    
    γB1_t = γB1_tilde(nf)
    γB2_t = γB2_tilde(nf)
    γB3_t = γB3_tilde(nf)  

    order1 = αs/(4π) * γB0_t
    order2 = (αs/(4π))^2 * γB1_t
    order3 = (αs/(4π))^3 * γB2_t
    order4 = (αs/(4π))^4 * γB3_t

    if order == 1 
        total = order1
    elseif order == 2  
        total = order1 + order2
    elseif order == 3 
        total = order1 + order2 + order3
    elseif order == 4 
        total = order1 + order2 + order3 + order4
    end

    return total
end

# numerical integration
#@everywhere function γJ_final(; μ::Float64, νdQ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
#
#    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
#    Γ = Γ_func(αs=αs, order=order+1, nf=nf)
#    γJ = γB_func(αs=αs, order=order, nf=nf) + γS_func(αs=αs, order=order, nf=nf)
#
#    γJ_f = 2*Γ*log(νdQ) + γJ
#
#    return γJ_f
#end