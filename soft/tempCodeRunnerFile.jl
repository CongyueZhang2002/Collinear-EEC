# TMD soft function up to 3 loops from https://arxiv.org/pdf/1604.01404 

using Distributed
include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")
include("..\\anomalous dim\\cusp\\cusp.jl")

# constants
@everywhere function cT1_func(nf)
    #cT1 = -2*CF*z2
    return -4.38648 #cT1
end

@everywhere function cT2_func(nf)
    #cT2 = (
    #    CA*CF*(-67*z2/3 - 154*z3/9 + 10*z4 + 2428/81)
    #    + CF*nf*(10*z2/3 + 28*z3/9 - 328/81)
    #)
    return -31.5376 #cT2
end

@everywhere function cT3_func(nf)
    #cT3 = (
    #    CF*CA^2 * (928/9*z3^2 + 1100/9*z2*z3 - 151132/243*z3 - 297481/729*z2 + 3649/27*z4 + 1804/9*z5 - 3086/27*z6 + 5211949/13122) 
    #    + CF*CA*nf * (40/9*z3*z2 + 74530/729*z2 + 8152/81*z3 - 416/27*z4 - 184/3*z5 - 412765/6561) 
    #    + CF^2*nf * (-80/3*z3*z2 + 275/9*z2 + 3488/81*z3 + 152/9*z4 + 224/9*z5 - 42727/486) 
    #    + CF*nf^2 * (-136/27*z2 - 560/243*z3 - 44/27*z4 - 256/6561)
    #)
    return -2002.91 #cT3
end

@everywhere function cs1_func(nf)
    #cs1 = 2*CF*z2
    return 4.38648 #cs1
end

@everywhere function cs2_func(nf)
    #cs2 = (
    #    CA*CF*(67*z2/9 - 22*z3/9 - 30*z4 + 2428/81)
    #    + CF*nf*(-10*z2/9 + 4*z3/9 - 328/81)
    #)
    return -8.36717 #cs2
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

    return 76.2110 #γs1
end

@everywhere function γs2_thre_func(nf)

    #γs2 = (
    #      CA^2*CF * (-176/3*z3*z2 + 12650/81*z2 + 1316/3*z3 - 176*z4 - 192*z5 - 136781/729) 
    #    + CA*CF*nf * (-2828/81*z2 - 728/27*z3 + 48*z4 + 11842/729) 
    #    + CF*CF*nf * (-4*z2 - 304/9*z3 - 16*z4 + 1711/27) 
    #    + CF*nf^2 * (40/27*z2 - 112/27*z3 + 2080/729)
    #)

    return 659.025 #γs2
end

# rapidity anomalous dimensions
@everywhere function γr0_func(nf)

    γr0 = 0

    return γr0
end

@everywhere function γr1_func(nf)

    #γr1 = CF*CA * (28*z3 - 808/27) + 112*CF*nf/27

    return 42.5813 #γr1
end

@everywhere function γr2_func(nf)
    
    γr2 = (
        CF*CA^2 * (-176/3*z3*z2 + 6392*z2/81 + 12328*z3/27 + 154*z4/3 - 192*z5 - 297029/729)
        + CF*CA * nf * (-824*z2/81 - 904*z3/27 + 20*z4/3 + 62626/729)
        + CF*nf^2 * (-32*z3/9 - 1856/729)
        + CF*CF*nf * (-304*z3/9 - 16*z4 + 1711/27)
    )

    return γr2
end
print(γr2_func(5))