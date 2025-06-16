using PolyLog
using Interpolations
using CSV 
using DataFrames
include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")

function A_func(x)

    setprecision(96)
    z=BigFloat(x)

    value = (
        (CF*(-3 + 2*z)*(3*(2 - 3*z)*z + 2*(3 - 6*z + 2*z^2)*log(1-z)))/(4*(-1 + z)*z^5)
    )

    return value
end

function B_func(x) 

    #setprecision(96)
    #z=BigFloat(x)
    z=x

    value = (
        (4*π^2*(-1266 + 3143*z - 2647*z^2 + 585*z^3 + 130*z^4 - 120*z^5 + 120*z^6) + 
    2*z*(9320 - 27552*z + 14966*z^2 + 902*z^3 - 17359*z^4 + 75748*z^5 - 115200*z^6 + 57600*z^7) + 
    4*(-4880 + 12412*z - 11322*z^2 + 3571*z^3 - 3225*z^4 + 31035*z^5 - 147846*z^6 + 321680*z^7 - 316800*z^8 + 115200*z^9)*log(1 - z) - 
    2*z*(11424 - 25029*z + 10971*z^2 - 742*z^3 + 18696*z^4 - 138600*z^5 + 412960*z^6 - 518400*z^7 + 230400*z^8)*log(z) - 
    3*(-1 + z)*z^(3/2)*(1435 + 547*z + 992*z^2 - 160*z^3 + 960*z^4)*(log((1 - sqrt(z))/(1 + sqrt(z)))*log(z) - 2*reli(2, -sqrt(z)) + 2*reli(2, sqrt(z))) + 
    24*(-1 + z)*(-952 + 1431*z - 315*z^2 + 40*z^3 + 340*z^4 - 2660*z^5 + 14680*z^6 - 28800*z^7 + 19200*z^8)*(reli(2, 1 - z) - reli(2, z)) - 
    4*(-314 + 760*z - 721*z^2 + 140*z^3 - 15*z^4 + 184*z^5 - 235*z^6 + 91*z^7)*(π^2 + 3*log(1 - z)^2 + 6*reli(2, z)) + 
    120*(-1 + z)*z*(3 + z^2 + 2*z^3 - z^4 + 2*z^5)*(3*log(-1 + z^(-1))*log((1 + sqrt(z))/(1 - sqrt(z)))^2 + 2*π^2*log(1 - z) - 
        24*reli(3, sqrt(z)/(-1 + sqrt(z))) - 24*reli(3, sqrt(z)/(1 + sqrt(z))) + 6*reli(3, z/(-1 + z))) + 
    120*(1 - 9*z + 9*z^2 - z^3 - z^4 + 3*z^5 - 3*z^6 + 2*z^7)*(log(1 - z)^3 + 6*log(1 - z)*reli(2, z) - 12*(reli(3, z) + reli(3, z/(-1 + z)))) + 
    15*(2*π^2*(7 + 71*z - 66*z^2 + 8*z^3) - 5*z*(2050 - 4115*z + 1825*z^2 + 48*z^3 - 1568*z^4 + 8852*z^5 - 14400*z^6 + 7200*z^7) - 
        2*(1801 - 4801*z + 3269*z^2 - 489*z^3 - 100*z^4 + 10960*z^5 - 77700*z^6 + 193040*z^7 - 198000*z^8 + 72000*z^9)*log(1 - z) + 
        4*z*(561 - 939*z + 428*z^2 + 10*z^3 + 1190*z^4 - 16650*z^5 + 60520*z^6 - 81000*z^7 + 36000*z^8)*log(z) + 
        15*(-1 + z)*z^(3/2)*(-1 + 3*z)*(log((1 - sqrt(z))/(1 + sqrt(z)))*log(z) - 2*reli(2, -sqrt(z)) + 2*reli(2, sqrt(z))) - 
        12*(-1 + z)*(-187 + 222*z - 72*z^2 - 920*z^5 + 7840*z^6 - 18000*z^7 + 12000*z^8)*(reli(2, 1 - z) - reli(2, z)) - 
        40*(-9 + 24*z - 18*z^2 + 4*z^3 + z^7)*(π^2 + 3*log(1 - z)^2 + 6*reli(2, z)) - 1440*(1 - z)*z^5*(1 - 16*z + 66*z^2 - 100*z^3 + 50*z^4)*
        ((-log(1 - z) + log(z))*(π^2/3 + log(1 - z)^2 + 2*reli(2, z)) + 6*(reli(3, z/(-1 + z)) - z3))) + 
    4*(-2*π^2*(4193 - 10159*z + 8812*z^2 - 2246*z^3 + 160*z^4 + 60*z^5 + 120*z^6) + 
        z*(63298 - 143577*z + 72305*z^2 + 2064*z^3 - 31000*z^4 + 157060*z^5 - 244800*z^6 + 122400*z^7) + 
        2*(-3007 + 9329*z - 11309*z^2 + 6201*z^3 - 2716*z^4 + 48122*z^5 - 283140*z^6 + 667280*z^7 - 673200*z^8 + 244800*z^9)*log(1 - z) - 
        2*z*(19938 - 38295*z + 17261*z^2 - 336*z^3 + 13052*z^4 - 126900*z^5 + 422480*z^6 - 550800*z^7 + 244800*z^8)*log(z) - 
        30*(-1 + z)*z^(3/2)*(-1 + 11*z)*(log((1 - sqrt(z))/(1 + sqrt(z)))*log(z) - 2*reli(2, -sqrt(z)) + 2*reli(2, sqrt(z))) + 
        12*(-1 + z)*(-3323 + 4726*z - 1126*z^2 + 160*z^3 + 320*z^4 - 4040*z^5 + 28480*z^6 - 61200*z^7 + 40800*z^8)*(reli(2, 1 - z) - reli(2, z)) + 
        20*(87 - 211*z + 296*z^2 - 96*z^3 + 25*z^4 - 17*z^5 + 10*z^6 + 4*z^7)*(π^2 + 3*log(1 - z)^2 + 6*reli(2, z)) + 
        120*z^5*(1 + z^2)*(log(1 - z)^3 - log(1 - z)*(π^2 - 6*reli(2, z)) - 12*reli(3, z)) + 
        240*(5 - 21*z + 18*z^2 - 4*z^3)*(log(1 - z)^3 + 6*log(1 - z)*reli(2, z) - 12*(reli(3, z) + reli(3, z/(-1 + z)))) + 
        2880*(1 - z)*z^5*(3 - 31*z + 116*z^2 - 170*z^3 + 85*z^4)*((-log(1 - z) + log(z))*(π^2/3 + log(1 - z)^2 + 2*reli(2, z)) + 
        6*(reli(3, z/(-1 + z)) - z3))) - 360*z^4*(3 - 42*z + 318*z^2 - 1196*z^3 + 2196*z^4 - 1920*z^5 + 640*z^6)*
    ((-log(1 - z) + log(z))*(π^2/3 + log(1 - z)^2 + 2*reli(2, z)) + 6*(reli(3, z/(-1 + z)) - z3)) - 
    360*z^4*(-1 + 2*z)*(1 - z + z^2)*(π^2*log(z) - 2*reli(3, z/(-1 + z)) - 16*z3))/(3240*(1 - z)*z^5)
        )

    return value
end

function Bzto0_func(x)

    setprecision(96)
    z=BigFloat(x)

    value = (
        -97340143/11907000 + 430013/(24300*z) - (522659869*z)/2976750 + (20076584289173*z^2)/8644482000 - (689192297878237073*z^3)/122717066472000 + 
        (918677086997*z^4)/205821000 + (10676*z2)/135 - (53*z2)/(27*z) - (67313*z*z2)/135 + (23111*z^2*z2)/27 + (292813*z^3*z2)/1890 - 
        (143204*z^4*z2)/135 + log(z)*(-186253/3780 - 98/(45*z) + z^3*(836161740053/340540200 - (13444*z2)/9) + z*(5233933/11340 - (2536*z2)/9) + 
        (86*z2)/3 + z^4*(-574991072/467775 + (6716*z2)/9) + z^2*(-87985166/51975 + 1028*z2) + 
        nf*(86501/18900 + 53/(360*z) + z^2*(15059174/51975 - 176*z2) + z^4*(64740629/294840 - (400*z2)/3) - (8*z2)/3 + 
            z*(-3967303/56700 + (128*z2)/3) + z^3*(-189580613/432432 + (800*z2)/3))) - (778*z3)/9 + (2*z3)/(9*z) + (7610*z*z3)/9 - 
        (27758*z^2*z3)/9 + (40330*z^3*z3)/9 - (20150*z^4*z3)/9 + 
        nf*(2987627/496125 - 4913/(5400*z) - (92*z2)/9 + z^3*(129450847890323/129859329600 - (202*z2)/9 - 800*z3) + 
        z*(970494187/71442000 + (760*z2)/9 - 128*z3) + 8*z3 + z^4*(-96409889267/120736980 + 192*z2 + 400*z3) + 
        z^2*(-540097312157/1440747000 - (1430*z2)/9 + 528*z3))
    )
    return value
end

function C_func(x)

    χ_list = range(0.9, stop=179.1, length=100)
    χ = acos(1-2*x)
    EEC_NNLO = df_NNLO[:, 2]/(sin(χ))*x

    interpolator = interpolate((χ_list,), EEC_NNLO, Gridded(Linear()))

    value = interpolator(χ*180/π)/x

    return value
end

# https://arxiv.org/pdf/1708.04093
function perturbation_func(; z::Float64, Q::Float64)

    order = nloops_FO

    αs = alpha_s_func(μren_ratio*Q)
    LQ = log(μren_ratio^2)

    β0 = β0_func(nf)
    β1 = β1_func(nf)

    x = z
    A = A_func(x)
    if x > 0.01
        B = B_func(x)
    else
        B = Bzto0_func(x)
    end
    C = C_func(x)

    order1 = (αs/(2π))*A
    order2 = (αs/(2π))^2*(
          B 
        + 1/2*β0*LQ*A
        )
    order3 = (αs/(2π))^3*(
          C 
        + β0*LQ*B
        + (1/4*β1*LQ + 1/4*β0^2*LQ^2)*A
        )    

    if order == 0 
        total = order1
    elseif order == 1 
        total = order1 + order2
    elseif order == 2 
        total = order1 + order2 + order3
    end

    return 2*total

end

function perturbation_NNLO_sigma_χ(; χ::Float64, Q::Float64)

    αs = alpha_s_func(Q)

    χ_list = range(0.9, stop=179.1, length=100)
    EEC_LO = df_LO[:, 2]
    EEC_NLO = df_NLO[:, 2]
    EEC_NNLO = df_NNLO[:, 2]

    as = (αs/(2π))
    EEC_list = as*EEC_LO + as^2*EEC_NLO + as^3*EEC_NNLO

    interpolator = interpolate((χ_list,), EEC_list, Gridded(Linear()))
    
    EEC = interpolator(χ*180/π)

    return EEC
end

function perturbation_sigma_χ(; χ::Float64, Q::Float64)

    order = nloops_FO

    if order < 2
        z = 0.5*(1-cos(χ))
        part = perturbation_func(z=z, Q=Q)
        total = 0.5*sin(χ)*part
    else
        total = perturbation_NNLO_sigma_χ(χ=χ, Q=Q)
    end
    
    return total
end