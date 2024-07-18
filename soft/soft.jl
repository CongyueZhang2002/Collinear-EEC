# TMD soft function up to 3 loops from https://arxiv.org/pdf/1604.01404 

using Distributed
include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")
include("..\\anomalous dim\\cusp\\cusp.jl")

# constants
@everywhere function cT1_func(nf)
    #cT1 = -2*CF*z2
    #return cT1
    return -4.38648
end

@everywhere function cT2_func(nf)
    #cT2 = (
    #    CA*CF*(-67*z2/3 - 154*z3/9 + 10*z4 + 2428/81)
    #    + CF*nf*(10*z2/3 + 28*z3/9 - 328/81)
    #)
    #return cT2
    return -31.5376
end

@everywhere function cT3_func(nf)
    #cT3 = (
    #    CF*CA^2 * (928/9*z3^2 + 1100/9*z2*z3 - 151132/243*z3 - 297481/729*z2 + 3649/27*z4 + 1804/9*z5 - 3086/27*z6 + 5211949/13122) 
    #    + CF*CA*nf * (40/9*z3*z2 + 74530/729*z2 + 8152/81*z3 - 416/27*z4 - 184/3*z5 - 412765/6561) 
    #    + CF^2*nf * (-80/3*z3*z2 + 275/9*z2 + 3488/81*z3 + 152/9*z4 + 224/9*z5 - 42727/486) 
    #    + CF*nf^2 * (-136/27*z2 - 560/243*z3 - 44/27*z4 - 256/6561)
    #)
    #return cT3
    return -2002.91
end

@everywhere function cs1_func(nf)
    #cs1 = 2*CF*z2
    #return cs1
    return 4.38648
end

@everywhere function cs2_func(nf)
    #cs2 = (
    #    CA*CF*(67*z2/9 - 22*z3/9 - 30*z4 + 2428/81)
    #    + CF*nf*(-10*z2/9 + 4*z3/9 - 328/81)
    #)
    #return cs2
    return -8.36717
end

# threshold soft anomalous dimensions
@everywhere function γs0_thre_func(nf)

    γs0 = 0

    return γs0
end

@everywhere function γs1_thre_func(nf)

    #γs1 = (
    #      CA*CF * (22/3*z2 + 28*z3 - 808/27) 
    #    + CF*nf * (112/27 - 4/3*z2)
    #)

    #return γs1
    return 76.2110
end

@everywhere function γs2_thre_func(nf)

    #γs2 = (
    #      CA^2*CF * (-176/3*z3*z2 + 12650/81*z2 + 1316/3*z3 - 176*z4 - 192*z5 - 136781/729) 
    #    + CA*CF*nf * (-2828/81*z2 - 728/27*z3 + 48*z4 + 11842/729) 
    #    + CF*CF*nf * (-4*z2 - 304/9*z3 - 16*z4 + 1711/27) 
    #    + CF*nf^2 * (40/27*z2 - 112/27*z3 + 2080/729)
    #)

    #return γs2
    return 659.025
end

# rapidity anomalous dimensions
@everywhere function γr0_func(nf)

    γr0 = 0

    return γr0
end

@everywhere function γr1_func(nf)

    #γr1 = CF*CA * (28*z3 - 808/27) + 112*CF*nf/27

    #return γr1
    return 42.5813
end

@everywhere function γr2_func(nf)
    
    #γr2 = (
    #    CF*CA^2 * (-176/3*z3*z2 + 6392*z2/81 + 12328*z3/27 + 154*z4/3 - 192*z5 - 297029/729)
    #    + CF*CA * nf * (-824*z2/81 - 904*z3/27 + 20*z4/3 + 62626/729)
    #    + CF*nf^2 * (-32*z3/9 - 1856/729)
    #    + CF*CF*nf * (-304*z3/9 - 16*z4 + 1711/27)
    #)

    #return γr2
    return 684.094
end

# Soft function (eqn 9S)
@everywhere function S1_func(Lb,Lr,nf)

    cT1 = cT1_func(nf)

    Γ0 = Γ0_func(nf)

    γs0 = γs0_thre_func(nf)
    γr0 = γr0_func(nf)

    S1 = cT1 + (1/2)*Γ0*Lb^2 + γr0*Lr - Lb*(γs0 + Γ0*Lr)

    return S1
end

@everywhere function S2_func(Lb,Lr,nf)

    cT1 = cT1_func(nf)
    cT2 = cT2_func(nf)

    β0 = β0_func(nf)  
    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)

    γs0 = γs0_thre_func(nf)
    γs1 = γs1_thre_func(nf)
    γr0 = γr0_func(nf)
    γr1 = γr1_func(nf)

    S2 = (
          cT2 
        + γr1*Lr 
        + (1/6)*Γ0*Lb^3*β0
        + Lb^2*(1/2*Γ1 - 1/2*γs0*β0 - 1/2*Γ0*Lr*β0) 
        + Lb*(-γs1 + cT1*β0 + Lr*(-Γ1 + γr0*β0))
    )

    return S2
end

@everywhere function S3_func(Lb,Lr,nf)

    cT1 = cT1_func(nf)
    cT2 = cT2_func(nf)
    cT3 = cT3_func(nf)

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)  

    γs0 = γs0_thre_func(nf)
    γs1 = γs1_thre_func(nf)
    γs2 = γs2_thre_func(nf)
    γr0 = γr0_func(nf)
    γr1 = γr1_func(nf)
    γr2 = γr2_func(nf)

    S3 = (
          cT3 
        + γr2*Lr  
        + (1/12)*Γ0*Lb^4*β0^2
        + Lb^3*(1/3*Γ1*β0 - 1/3*γs0*β0^2 - 1/3*Γ0*Lr*β0^2 + Γ0*β1/6) 
        + Lb^2*(1/2*Γ2 - γs1*β0 + cT1*β0^2 - γs0*β1/2 + Lr*(-Γ1*β0 + γr0*β0^2 - Γ0*β1/2)) 
        + Lb*(-γs2 + 2*cT2*β0 + cT1*β1 + Lr*(-Γ2 + 2*γr1*β0 + γr0*β1)) 
    )

    return S3
end

@everywhere function S_func(; b::Float64, μS::Float64, νS::Float64, αs::Float64, order::Int64, nf::Int64)

    Lb = log((b*μS/b0)^2)
    Lr = log((b*νS/b0)^2)

    S1 = S1_func(Lb,Lr,nf)
    S2 = S2_func(Lb,Lr,nf)
    S3 = S3_func(Lb,Lr,nf)

    if order == 1 
        total = (αs/4π)*S1
    elseif order == 2  
        total = (αs/4π)*S1 + (αs/4π)^2*S2
    elseif order == 3 
        total = (αs/4π)*S1 + (αs/4π)^2*S2 + (αs/4π)^3*S3
    end

    return exp(total)
end