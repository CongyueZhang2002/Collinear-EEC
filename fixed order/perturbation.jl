# https://inspirehep.net/files/57bd06a70eeacbc1cf3f5fd908b9318b 
using PolyLog
include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")

#------------------------------------------

function g11_func(z)
    value = log(1-z)
    return value
end

function g12_func(z)
    value = log(z)
    return value
end

#------------------------------------------

function g21_func(z)
    value = log(1 - z)^2 + 2*(reli(2, z) + z2)
    return value
end

function g22_func(z)
    value = reli(2, 1 - z) - reli(2, z)
    return value
end

function g23_func(z)
    value = (
        log((1 - sqrt(z))/(1 + sqrt(z)))*log(z) - 
        2*reli(2, -sqrt(z)) + 2*reli(2, sqrt(z))
    )
    return value
end

function g24_func(z)
    value = z2
    return value
end

#function g25_func(z)
#    value = (
#    (-im)*(-(log(sqrt(z/(1 - z)))*log((1 + im*sqrt(z/(1 - z)))/
#    (1 - im*sqrt(z/(1 - z))))) - reli(2, (-im)*sqrt(z/(1 - z))) + 
#    reli(2, im*sqrt(z/(1 - z))))
#    )
#    return value
#end

#------------------------------------------

function g31_func(z)
    value = (
        -((-log(1 - z) + log(z))*(log(1 - z)^2 + 
        2*(reli(2, z) + z2))) - 6*(reli(3, -(z/(1 - z))) - 
       z3)
    )
    return value
end

function g32_func(z)
    value = (
        log(1 - z)^3 + 6*log(1 - z)*reli(2, z) - 
        12*(reli(3, z) + reli(3, -(z/(1 - z))))
    )
    return value
end

function g33_func(z)
    value = (
        log(1 - z)^3 - 12*reli(3, z) + 
        6*log(1 - z)*(reli(2, z) - z2)
    )
    return value
end

function g34_func(z)
    value = (
        reli(3, -(z/(1 - z))) - 3*log(z)*z2 + 8*z3
    )
    return value
end

function g35_func(z)
    value = (
        log((1 + sqrt(z))/(1 - sqrt(z)))^2*log((1 - z)/z) - 
       8*(reli(3, -(sqrt(z)/(1 - sqrt(z)))) + reli(3, 
          sqrt(z)/(1 + sqrt(z)))) + 2*reli(3, -(z/(1 - z))) + 
       4*log(1 - z)*z2
    )
    return value
end

function g36_func(z)
    value = (
        log(1 - z)^3 - 15*log(1 - z)*z2
    )
    return value
end

function g37_func(z)
    value = (
        log(1 - z)*(log(1 - z)*log(z) + reli(2, z) - 
        (15*z2)/2)
    )
    return value
end

function g38_func(z)
    value = z3
    return value
end

#function g39_func(z)
#    value = (
#        -(π^2*(log((1 - (im*sqrt(z))/sqrt(1 - z))/2) + 
#          log((1 + (im*sqrt(z))/sqrt(1 - z))/2))) + 
#       2*(log((1 - (im*sqrt(z))/sqrt(1 - z))/2)^3 + 
#         log((1 + (im*sqrt(z))/sqrt(1 - z))/2)^3) - 
#       3*(-log(1 - (im*sqrt(z))/sqrt(1 - z)) + 
#          log(1 + (im*sqrt(z))/sqrt(1 - z)))^2*((-im)*π + 
#         2*log((im*sqrt(z))/sqrt(1 - z))) + 3*reli(3, -(z/(1 - z))) - 
#       12*(-reli(3, (1 - (im*sqrt(z))/sqrt(1 - z))/2) - 
#         reli(3, (1 + (im*sqrt(z))/sqrt(1 - z))/2) + 
#         reli(3, ((-im)*sqrt(z))/sqrt(1 - z)) + reli(3, 
#          (im*sqrt(z))/sqrt(1 - z)) + reli(3, (-2*sqrt(z))/
#           ((im - sqrt(z)/sqrt(1 - z))*sqrt(1 - z))) + 
#         reli(3, (2*sqrt(z))/((im + sqrt(z)/sqrt(1 - z))*sqrt(1 - z))) - 
#         z3)
#    )
#    return value
#end

#------------------------------------------

function A1_func(z)
    value = (CF*(-3 + 2*z)*(3*(2 - 3*z)*z + 2*(3 - 6*z + 2*z^2)*log(1 - z)))/(4*(-1 + z)*z^5)
    return value
end

function Blc_func(z)

    g11 = g11_func(z)
    g12 = g12_func(z)
    g21 = g21_func(z)
    g22 = g22_func(z)
    g23 = g23_func(z)
    g24 = g24_func(z)

    g31 = g31_func(z)
    g32 = g32_func(z)
    g33 = g33_func(z)
    g34 = g34_func(z)
    g35 = g35_func(z)
    g36 = g36_func(z)
    g37 = g37_func(z)
    g38 = g38_func(z)


    value = ( 
        (63298 - 143577*z + 72305*z^2 + 2064*z^3 - 31000*z^4 + 157060*z^5 - 
        244800*z^6 + 122400*z^7)/(1440*(1 - z)*z^4) - 
      ((3007 - 9329*z + 11309*z^2 - 6201*z^3 + 2716*z^4 - 48122*z^5 + 
         283140*z^6 - 667280*z^7 + 673200*z^8 - 244800*z^9)*g11)/
       (720*(1 - z)*z^5) - ((19938 - 38295*z + 17261*z^2 - 336*z^3 + 
         13052*z^4 - 126900*z^5 + 422480*z^6 - 550800*z^7 + 244800*z^8)*
        g12)/(720*(1 - z)*z^4) + 
      ((87 - 211*z + 296*z^2 - 96*z^3 + 25*z^4 - 17*z^5 + 10*z^6 + 4*z^7)*
        g21)/(24*(1 - z)*z^5) + 
      ((3323 - 4726*z + 1126*z^2 - 160*z^3 - 320*z^4 + 4040*z^5 - 28480*z^6 + 
         61200*z^7 - 40800*z^8)*g22)/(120*z^5) - 
      ((1 - 11*z)*g23)/(48*z^(7/2)) - 
      ((4193 - 10159*z + 8812*z^2 - 2246*z^3 + 160*z^4 + 60*z^5 + 120*z^6)*
        g24)/(120*(1 - z)*z^5) - 2*(3 - 31*z + 116*z^2 - 170*z^3 + 85*z^4)*
       g31 + ((5 - 21*z + 18*z^2 - 4*z^3)*g32)/(6*(1 - z)*z^5) + 
      ((1 + z^2)*g33)/(12*(1 - z))
    )
    return value
end

function Bnlc_func(z)

    g11 = g11_func(z)
    g12 = g12_func(z)
    g21 = g21_func(z)
    g22 = g22_func(z)
    g23 = g23_func(z)
    g24 = g24_func(z)

    g31 = g31_func(z)
    g32 = g32_func(z)
    g33 = g33_func(z)
    g34 = g34_func(z)
    g35 = g35_func(z)
    g36 = g36_func(z)
    g37 = g37_func(z)
    g38 = g38_func(z)


    value = ( 
        (9320 - 27552*z + 14966*z^2 + 902*z^3 - 17359*z^4 + 75748*z^5 - 
        115200*z^6 + 57600*z^7)/(720*(1 - z)*z^4) - 
      ((4880 - 12412*z + 11322*z^2 - 3571*z^3 + 3225*z^4 - 31035*z^5 + 
         147846*z^6 - 321680*z^7 + 316800*z^8 - 115200*z^9)*g11)/
       (360*(1 - z)*z^5) - ((11424 - 25029*z + 10971*z^2 - 742*z^3 + 
         18696*z^4 - 138600*z^5 + 412960*z^6 - 518400*z^7 + 230400*z^8)*
        g12)/(720*(1 - z)*z^4) + 
      ((314 - 760*z + 721*z^2 - 140*z^3 + 15*z^4 - 184*z^5 + 235*z^6 - 91*z^7)*
        g21)/(120*(1 - z)*z^5) + 
      ((952 - 1431*z + 315*z^2 - 40*z^3 - 340*z^4 + 2660*z^5 - 14680*z^6 + 
         28800*z^7 - 19200*z^8)*g22)/(60*z^5) + 
      ((1435 + 547*z + 992*z^2 - 160*z^3 + 960*z^4)*g23)/(480*z^(7/2)) - 
      ((1266 - 3143*z + 2647*z^2 - 585*z^3 - 130*z^4 + 120*z^5 - 120*z^6)*
        g24)/(60*(1 - z)*z^5) + 
      ((3 - 42*z + 318*z^2 - 1196*z^3 + 2196*z^4 - 1920*z^5 + 640*z^6)*
        g31)/(4*(1 - z)*z) + ((1 - 9*z + 9*z^2 - z^3 - z^4 + 3*z^5 - 
         3*z^6 + 2*z^7)*g32)/(12*(1 - z)*z^5) - 
      ((1 - 2*z)*(1 - z + z^2)*g34)/(2*(1 - z)*z) - 
      ((3 + z^2 + 2*z^3 - z^4 + 2*z^5)*g35)/(4*z^4)
    )
    return value
end    

print(Bnlc_func(0.5))

function BNf_func(z)

    g11 = g11_func(z)
    g12 = g12_func(z)
    g21 = g21_func(z)
    g22 = g22_func(z)
    g23 = g23_func(z)
    g24 = g24_func(z)

    g31 = g31_func(z)
    g32 = g32_func(z)
    g33 = g33_func(z)
    g34 = g34_func(z)
    g35 = g35_func(z)
    g36 = g36_func(z)
    g37 = g37_func(z)
    g38 = g38_func(z)


    value = ( 
        -(2050 - 4115*z + 1825*z^2 + 48*z^3 - 1568*z^4 + 8852*z^5 - 14400*z^6 + 
        7200*z^7)/(144*(1 - z)*z^4) - 
     ((1801 - 4801*z + 3269*z^2 - 489*z^3 - 100*z^4 + 10960*z^5 - 77700*z^6 + 
        193040*z^7 - 198000*z^8 + 72000*z^9)*g11)/(360*(1 - z)*z^5) + 
     ((561 - 939*z + 428*z^2 + 10*z^3 + 1190*z^4 - 16650*z^5 + 60520*z^6 - 
        81000*z^7 + 36000*z^8)*g12)/(180*(1 - z)*z^4) + 
     ((9 - 24*z + 18*z^2 - 4*z^3 - z^7)*g21)/(6*(1 - z)*z^5) - 
     ((187 - 222*z + 72*z^2 + 920*z^5 - 7840*z^6 + 18000*z^7 - 12000*z^8)*
       g22)/(60*z^5) + ((1 - 3*z)*g23)/(48*z^(7/2)) + 
     ((7 + 71*z - 66*z^2 + 8*z^3)*g24)/(60*(1 - z)*z^5) + 
     2*(1 - 16*z + 66*z^2 - 100*z^3 + 50*z^4)*g31
    )
    return value
end

function Bqq_func(z)

    g11 = g11_func(z)
    g12 = g12_func(z)
    g21 = g21_func(z)
    g22 = g22_func(z)
    g23 = g23_func(z)
    g24 = g24_func(z)

    g31 = g31_func(z)
    g32 = g32_func(z)
    g33 = g33_func(z)
    g34 = g34_func(z)
    g35 = g35_func(z)
    g36 = g36_func(z)
    g37 = g37_func(z)
    g38 = g38_func(z)


    value = ( 
        -(13696 - 7735*z - 3001*z^2 - 1028*z^3 + 24394*z^4 - 104356*z^5 + 
        158400*z^6 - 79200*z^7)/(1440*(1 - z)*z^4) + 
     ((1637 - 2812*z + 1146*z^2 + 28*z^3 - 2880*z^4 + 20096*z^5 - 83680*z^6 + 
        138600*z^7 - 79200*z^8)*g11)/(360*z^5) + 
     ((5358 - 4239*z + 1448*z^2 + 231*z^3 - 7253*z^4 + 48870*z^5 - 
        143080*z^6 + 178200*z^7 - 79200*z^8)*g12)/(360*(1 - z)*z^4) - 
     ((325 - 748*z + 520*z^2 - 100*z^3 - 135*z^4 + 220*z^5 - 202*z^6)*
       g21)/(240*z^5) - ((113 + 225*z + 41*z^2 + 75*z^3 + 345*z^4 - 
        1930*z^5 + 10280*z^6 - 19800*z^7 + 13200*z^8)*g22)/(60*z^5) + 
     ((240 + 320*z + 375*z^2 + 496*z^3 - 80*z^4 + 480*z^5)*g23)/
      (240*z^(9/2)) + ((551 - 969*z + 840*z^2 - 592*z^3 + 260*z^4 - 240*z^5 + 
        240*z^6)*g24)/(120*(1 - z)*z^5) + ((96 - 68*z + 15*z^2)*g25)/
      (4*sqrt(1 - z)*z^(9/2)) + ((3 - 31*z + 199*z^2 - 636*z^3 + 880*z^4 - 
        440*z^5)*g31)/(4*z) - ((1 - 2*z)*(1 - z + z^2)*
       (g32/6 + g34 + g35/2))/(2*(1 - z)*z) + 
     ((22 - 11*z + z^2)*(-g32/2 + g35 - g36/2 + g37))/
      (8*z^5) + ((24 - 33*z + 12*z^2 - z^3)*(-3*g35 - 33*g38 + 
        g39))/(24*(1 - z)*z^5)
    )
    return value
end

#------------------------------------------

function perturbation_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    β0 = β0_func(nf)

    A = A1_func(z)
    Blc = Blc_func(z)
    Bnlc = Bnlc_func(z)
    BNf = BNf_func(z)

    B = (CF^2*Blc + CF*(CA-2CF)*Bnlc + CF*nf*TF*BNf)

    order1 = (αs/(2π))*A
    order2 = (αs/(2π))^2*B

    if order == 1 
        total = order1
    elseif order == 2 
        total = order1 + order2
    end

    return 2*real(total)

end

function perturbation_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end

#------------------------------------------

#function Blc_to1_func(z) 
#    value = (
#        -16129/288 + (53*log(1 - z)^2)/4 + 2*log(1 - z)^3 + 
#     (2893*z2)/40 + log(1 - z)*(-319/36 + 25*z2) + 
#     (1 - z)^2*(286777/288 + (239*log(1 - z)^2)/4 + (27*log(1 - z)^3)/2 + 
#       (356239*z2)/960 + log(1 - z)*(-47383/72 + 545*z2) - 
#       1230*z3) - 8*z3 + (-25/16 + (49*log(1 - z)^2)/12 + 
#       log(1 - z)^3/2 + (17*z2)/2 + log(1 - z)*(59/18 + 2*z2) + 
#       2*z3)/(1 - z) + (1 - z)*(-44665/288 + (175*log(1 - z)^2)/6 + 
#       (79*log(1 - z)^3)/12 + log(1 - z)*(757/3 - 85*z2) - 
#       (5197*z2)/30 + 449*z3) + 
#     (1 - z)^3*(-165863/72 + (449*log(1 - z)^2)/4 + (45*log(1 - z)^3)/2 + 
#       log(1 - z)*(60775/48 - 545*z2) - (46315*z2)/192 + 
#       2310*z3)
#    )
#    return value
#end

#function Bnlc_to1_func(z) 
#    value = (
#        -3305/72 + (27*log(1 - z)^2)/8 + log(1 - z)^3/2 + 
#        (6347*z2)/80 - 21*log(2)*z2 + 
#        log(1 - z)*(-2011/72 + 22*z2) + 
#        (1 - z)^2*(100897/108 + (493*log(1 - z)^2)/24 + (59*log(1 - z)^3)/12 + 
#          (54695*z2)/128 - 105*log(2)*z2 + 
#          log(1 - z)*(-18607/24 + (1007*z2)/2) - (5009*z3)/4) - 
#        (137*z3)/4 + (-35/16 + (11*log(1 - z)^2)/12 + 
#          log(1 - z)*(-35/72 + z2/2) + (11*z2)/4 + (3*z3)/2)/
#         (1 - z) + (1 - z)*(-7663/36 + (179*log(1 - z)^2)/24 + 2*log(1 - z)^3 + 
#          log(1 - z)*(671/3 - (237*z2)/2) - (8075*z2)/48 - 
#          42*log(2)*z2 + 464*z3) + 
#        (1 - z)^3*(-489985/192 + (731*log(1 - z)^2)/16 + (115*log(1 - z)^3)/12 + 
#          log(1 - z)*(68317/72 - 551*z2) + (14215*z2)/192 - 
#          198*log(2)*z2 + 2152*z3)
#    )
#    return value
#end

#function BNf_to1_func(z) 
#    value = (
#         2099/144 - log(1 - z)^2/2 + log(1 - z)*(361/36 - 4*z2) + 
#        (3/4 + log(1 - z)/18 - log(1 - z)^2/3 - z2)/(1 - z) - 
#        (1747*z2)/120 + (1 - z)^3*(62657/48 + (20*log(1 - z)^2)/3 + 
#          (14177*z2)/192 + log(1 - z)*(-11771/18 + 400*z2) - 
#          1200*z3) + (1 - z)*(1637/144 - (5*log(1 - z)^2)/3 + 
#          (7891*z2)/60 + log(1 - z)*(-1177/12 + 64*z2) - 
#          192*z3) + 12*z3 + (1 - z)^2*(-269485/432 - log(1 - z)^2/6 + 
#          log(1 - z)*(2647/6 - 264*z2) - (196571*z2)/960 + 792*z3)
#    )
#    return value
#end

#------------------------------------------

#function perturbation_to1_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)
#
#    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order, nf=nf)
#
#    A = A1_func(z)
#    Blc = Blc_to1_func(z)
#    Bnlc = Bnlc_to1_func(z)
#    BNf = BNf_to1_func(z)
#
#    B = CF^2*Blc + CF*(CA-2CF)*Bnlc + CF*nf*TF*BNf
#
#    order1 = (αs/(2π))*A
#    order2 = (αs/(2π))^2*B
#
#    if order == 1 
#        total = order1
#    end
#    if order == 2 
#        total = order1 + order2
#    end
#
#    return 2*total
#
#end

function perturbation_to1_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = perturbation_to1_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total

end