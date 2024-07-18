# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

using Distributed
using HCubature

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")
include("gammav.jl")

# Soft

@everywhere function γS0_func(nf)

    γS0 = 0.0

    return γS0
end

@everywhere function γS1_func(nf) 

    #β0 = β0_func(nf) 

    #γS1 = 2CF * (CA*(-64/9 + 28z3) + β0*(-56/9 + 2z2))

    #return γS1
    return 152.422
end

@everywhere function γS2_func(nf)

    #β0 = β0_func(nf)
    #β1 = β1_func(nf)

    #γS2 = 2CF * (
    #      CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
    #    + CA*β0 * (4697/54 + 484/27*z2 + 220/9*z3 - 112z4) 
    #    + β0^2 * (520/81 + 10/3*z2 - 28/3*z3) 
    #    + β1 * (-1711/54 + 2z2 + 152/9*z3 + 8z4)
    #)

    #return γS2
    return 1318.05
end

@everywhere function γS3_func(nf)

    #γν3 = γν3_func(nf)

    #γS3 = 4 * (
    #      CF^2*CF*nf * z2 
    #    + CF*nf^3 * (8/81*z2 - 880/243*z3 - 52/27*z4 + 256/6561) 
    #    + CA*CF*nf^2 * (56/9*z3*z2 + 5353/1458*z2 + 616/243*z3 - 26*z4 + 136/9*z5 + 166639/2187) 
    #    + CA*CF*CF*nf * (88/3*z3*z2 - 539/9*z2 + 16016/81*z3 + 1394/9*z4 + 1232/9*z5 - 528269/972) 
    #    + CA^3*CF * (7480/9*z3^2 + 242/9*z2*z3 - 486706/243*z3 - 1273/729*z2 - 128191/54*z4 - 2662/9*z5 + 48499/27*z6 + 66247055/26244) 
    #    + CF*CF*nf^2 * (-16/3*z3*z2 + 86/9*z2 - 2912/81*z3 - 152/9*z4 - 224/9*z5 + 46663/486) 
    #    + CA^2*CF*nf * (-1360/9*z3^2 - 352/9*z2*z3 + 111724/243*z3 - 48257/1458*z2 + 2017/3*z4 - 88/3*z5 - 8818/27*z6 - 648094/729)
    #) + γν3

    #return γS3
    return 10698.8
end

@everywhere function γS0_tilde(nf)
    return -γS0_func(nf)
end

@everywhere function γS1_tilde(nf)
    return -γS1_func(nf)
end

@everywhere function γS2_tilde(nf)
    return -γS2_func(nf)
end

@everywhere function γS3_tilde(nf)
    return -γS3_func(nf)
end

@everywhere function γS_tilde(; αs::Float64, order::Int64, nf::Int64)

    γS0_t = γS0_tilde(nf) 
    γS1_t = γS1_tilde(nf)
    γS2_t = γS2_tilde(nf)
    γS3_t = γS3_tilde(nf) 

    order1 = αs/(4π) * γS0_t
    order2 = (αs/(4π))^2 * γS1_t
    order3 = (αs/(4π))^3 * γS2_t
    order4 = (αs/(4π))^4 * γS3_t

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

#numerical integration
#@everywhere function γS_final(; μ::Float64, ν::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
#
#    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
#    Γ = Γ_func(αs=αs, order=order+1, nf=nf)
#    γS = γS_tilde(αs=αs, order=order, nf=nf)
#
#    γS_f = 4*Γ*log(μ/ν) + γS
#    
#    return γS_f
#end