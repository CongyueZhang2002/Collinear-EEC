#https://www.sciencedirect.com/science/article/pii/S0550321306004330
#timelike singlet splitting functions

using QuadGK
using PolyLog
include("harmonic polylogs.jl")
include("constants.jl")

#quark---------------------------------------------------------------------------------------------

function Hq1_Reg_func(x)

    value = 1/2*CF*(
        5 - 3*x - 2*(1 + x)*log(1-x) - 4*(1 + x)*log(x)
    )

    return value
end

function Hq1_Delta_func(x)

    value = 1/2*CF*(
        -9 + (4*π^2)/3
    )

    return value
end

function Hq1_D0_func(x)

    value = 1/2*CF*(
        -3 + 8*log(x)
    )

    return value
end

function Hq1_D1_func(x)

    value = 1/2*CF*(4)

    return value
end

function Hq2_Reg_func(x,nf)

    logx = log(x)    
    log1mx = log(1-x)
    log1px = log(1+x)

    value = (
        (16*(1 + x)*(11 + 19*x^2 + 6*nf*(-1 + x^2))*log1mx*reli(2, 1 - x))/(9*(-1 + x^2)) + 
        (2*(288*nf - 72*x + 864*nf*x - 216*x^2 + 576*nf*x^2 + 72*x^3 - 576*nf*x^3 + 72*x^4 - 864*nf*x^4 - 288*nf*x^5 + 144*x^6 - 144*x*log1mx + 144*x^2*log1mx - 144*x^3*log1mx + 144*x^4*log1mx - 288*x*logx - 
           432*x^3*logx - 144*x^4*logx)*reli(2, -x))/(81*x*(-1 + x^2)) + (2*(576*nf - 4068*x + 288*nf*x + 1980*x^2 - 1584*nf*x^2 + 2556*x^3 - 432*nf*x^3 - 3492*x^4 + 1152*nf*x^4 + 288*nf*x^5 + 720*x*log1mx + 
           720*x^2*log1mx + 720*x^3*log1mx + 720*x^4*log1mx + 3456*x*logx + 3456*x^2*logx + 3312*x^3*logx + 3312*x^4*logx)*reli(2, x))/(81*x*(-1 + x^2)) + 
        (2*(144*x - 144*x^2 + 144*x^3 - 144*x^4)*reli(3, (1 - x)/2))/(81*x*(-1 + x^2)) + (2*(-936*x + 432*nf*x - 648*x^2 + 432*nf*x^2 - 1512*x^3 - 432*nf*x^3 - 1224*x^4 - 432*nf*x^4)*reli(3, 1 - x))/(81*x*(-1 + x^2)) + 
        (2*(432*x + 144*x^2 + 864*x^3 + 288*x^4)*reli(3, -x))/(81*x*(-1 + x^2)) + (2*(-6192*x + 1296*nf*x - 6048*x^2 + 1296*nf*x^2 - 7776*x^3 - 1296*nf*x^3 - 7632*x^4 - 1296*nf*x^4)*reli(3, x))/(81*x*(-1 + x^2)) + 
        (2*(-144*x + 144*x^2 - 144*x^3 + 144*x^4)*reli(3, (2*x)/(-1 + x)))/(81*x*(-1 + x^2)) + (2*(144*x - 144*x^2 + 144*x^3 - 144*x^4)*reli(3, x/(1 + x)))/(81*x*(-1 + x^2)) + 
        (2*(-144*x + 144*x^2 - 144*x^3 + 144*x^4)*reli(3, (2*x)/(1 + x)))/(81*x*(-1 + x^2)) + (2*(144*x - 144*x^2 + 144*x^3 - 144*x^4)*reli(3, (1 + x)/2))/(81*x*(-1 + x^2)) + 
        (2*(592*nf - 6879*x + 3550*nf*x + 492*π^2*x - 24*nf*π^2*x + 429*x^2 - 3762*nf*x^2 + 312*π^2*x^2 + 264*nf*π^2*x^2 + 7023*x^3 - 4142*nf*x^3 - 240*π^2*x^3 - 429*x^4 + 3170*nf*x^4 - 84*π^2*x^4 - 288*nf*π^2*x^4 - 144*x^5 + 
           592*nf*x^5 + 24*π^2*x^6 - 48*x*log2^3 + 48*x^2*log2^3 - 48*x^3*log2^3 + 48*x^4*log2^3 + 12*π^2*x*log4 - 12*π^2*x^2*log4 + 12*π^2*x^3*log4 - 12*π^2*x^4*log4 - 144*nf*log1mx + 1962*x*log1mx + 
           1500*nf*x*log1mx - 132*π^2*x*log1mx + 1638*x^2*log1mx - 804*nf*x^2*log1mx - 156*π^2*x^2*log1mx - 1962*x^3*log1mx - 2220*nf*x^3*log1mx - 84*π^2*x^3*log1mx - 1638*x^4*log1mx + 
           948*nf*x^4*log1mx - 108*π^2*x^4*log1mx + 720*nf*x^5*log1mx + 72*x*log2^2*log1mx - 72*x^2*log2^2*log1mx + 72*x^3*log2^2*log1mx - 72*x^4*log2^2*log1mx - 144*nf*log1mx^2 - 
           1458*x*log1mx^2 - 72*nf*x*log1mx^2 - 882*x^2*log1mx^2 + 288*nf*x^2*log1mx^2 + 1458*x^3*log1mx^2 + 216*nf*x^3*log1mx^2 + 882*x^4*log1mx^2 - 144*nf*x^4*log1mx^2 - 144*nf*x^5*log1mx^2 - 
           24*x*log8*log1mx^2 + 24*x^2*log8*log1mx^2 - 24*x^3*log8*log1mx^2 + 24*x^4*log8*log1mx^2 + 312*x*log1mx^3 + 264*x^2*log1mx^3 - 264*x^3*log1mx^3 - 312*x^4*log1mx^3 - 288*nf*logx - 
           4086*x*logx + 4788*nf*x*logx + 84*π^2*x*logx - 144*nf*π^2*x*logx - 1710*x^2*logx + 3276*nf*x^2*logx + 108*π^2*x^2*logx - 144*nf*π^2*x^2*logx + 594*x^3*logx - 3588*nf*x^3*logx + 276*π^2*x^3*logx + 
           144*nf*π^2*x^3*logx - 1638*x^4*logx - 2412*nf*x^4*logx + 300*π^2*x^4*logx + 144*nf*π^2*x^4*logx - 144*x^5*logx - 624*nf*x^5*logx - 1440*x*log1mx*logx + 3744*x^2*log1mx*logx + 
           3168*x^3*log1mx*logx - 2016*x^4*log1mx*logx - 216*x*log1mx^2*logx - 216*nf*x*log1mx^2*logx - 72*x^2*log1mx^2*logx - 216*nf*x^2*log1mx^2*logx + 72*x^3*log1mx^2*logx + 
           216*nf*x^3*log1mx^2*logx + 216*x^4*log1mx^2*logx + 216*nf*x^4*log1mx^2*logx - 576*nf*logx^2 + 729*x*logx^2 - 288*nf*x*logx^2 - 1179*x^2*logx^2 + 936*nf*x^2*logx^2 - 2511*x^3*logx^2 + 
           108*nf*x^3*logx^2 - 531*x^4*logx^2 - 396*nf*x^4*logx^2 + 144*nf*x^5*logx^2 - 72*x^6*logx^2 + 108*x*log1mx*logx^2 + 108*x^2*log1mx*logx^2 + 36*x^3*log1mx*logx^2 + 36*x^4*log1mx*logx^2 + 
           264*x*logx^3 - 396*nf*x*logx^3 + 156*x^2*logx^3 - 396*nf*x^2*logx^3 + 792*x^3*logx^3 + 396*nf*x^3*logx^3 + 684*x^4*logx^3 + 396*nf*x^4*logx^3 - 72*π^2*x*log1px + 72*π^2*x^2*log1px - 
           144*π^2*x^3*log1px + 144*π^2*x^4*log1px + 72*x*log2^2*log1px - 72*x^2*log2^2*log1px + 72*x^3*log2^2*log1px - 72*x^4*log2^2*log1px + 288*nf*logx*log1px - 72*x*logx*log1px + 
           864*nf*x*logx*log1px - 216*x^2*logx*log1px + 576*nf*x^2*logx*log1px + 72*x^3*logx*log1px - 576*nf*x^3*logx*log1px + 72*x^4*logx*log1px - 864*nf*x^4*logx*log1px - 
           288*nf*x^5*logx*log1px + 144*x^6*logx*log1px - 144*x*logx^2*log1px + 144*x^2*logx^2*log1px - 72*x^3*logx^2*log1px + 72*x^4*logx^2*log1px - 12*x*log(64)*log1px^2 + 
           12*x^2*log(64)*log1px^2 - 12*x^3*log(64)*log1px^2 + 12*x^4*log(64)*log1px^2 - 144*x*(-1 + x - 3*x^2 + 3*x^3)*reli(3, 1 + x) + 9324*x*z3 - 1296*nf*x*z3 + 9468*x^2*z3 - 1296*nf*x^2*z3 + 
           4860*x^3*z3 + 1296*nf*x^3*z3 + 5292*x^4*z3 + 1296*nf*x^4*z3))/(81*x*(-1 + x^2))       
    )

    return value
end

function Hq2_Delta_func(x,nf)

    value = (
        (-62115 + 4570*nf + 9780*π^2 - 760*nf*π^2 + 106*π^4 + 480*(27 + nf)*z3)/270
    )

    return value
end

function Hq2_D0_func(x,nf)

    value = (
        (-2*(7629 - 396*π^2 + nf*(-494 + 24*π^2) - 5904*z3))/81
    )

    return value
end

function Hq2_D1_func(x,nf)

    value = (
        (-4*(-777 + 58*nf + 4*π^2))/27
    )

    return value
end

function Hq2_D2_func(x,nf)

    value = (
        (8*(-69 + 2*nf))/9
    )

    return value
end

function Hq2_D3_func(x,nf)

    value = (
        128/9
    )

    return value
end

function Hq_integrand(; x::Float64, α::Float64, order::Int64, fx::Float64, f1::Float64, nf::Int64)

    integrand0 = 1/2*f1

    Hq_Reg = Hq1_Reg_func(x)
    Hq_Delta = Hq1_Delta_func(x)
    Hq_D0_x = Hq1_D0_func(x)
    Hq_D0_1 = Hq1_D0_func(1)
    Hq_D1_x = Hq1_D1_func(x)
    Hq_D1_1 = Hq1_D1_func(1)

    integrand1 = (
          Hq_Reg*fx
        + Hq_Delta*f1 
        + (Hq_D0_x*fx - Hq_D0_1*f1)/(1-x)
        + (Hq_D1_x*fx - Hq_D1_1*f1)*log(1-x)/(1-x)
    )

    Hq_Reg = Hq2_Reg_func(x,nf)
    Hq_Delta = Hq2_Delta_func(x,nf)
    Hq_D0_x = Hq2_D0_func(x,nf)
    Hq_D0_1 = Hq2_D0_func(1,nf)
    Hq_D1_x = Hq2_D1_func(x,nf)
    Hq_D1_1 = Hq2_D1_func(1,nf)
    Hq_D2_x = Hq2_D2_func(x,nf)
    Hq_D2_1 = Hq2_D2_func(1,nf)
    Hq_D3_x = Hq2_D3_func(x,nf)
    Hq_D3_1 = Hq2_D3_func(1,nf)

    integrand2 = (
          Hq_Reg*fx
        + Hq_Delta*f1 
        + (Hq_D0_x*fx - Hq_D0_1*f1)/(1-x)
        + (Hq_D1_x*fx - Hq_D1_1*f1)*log(1-x)/(1-x)
        + (Hq_D2_x*fx - Hq_D2_1*f1)*(log(1-x))^2/(1-x)
        + (Hq_D3_x*fx - Hq_D3_1*f1)*(log(1-x))^3/(1-x)
    )/2

    if order < 0
        total = 0
    elseif order == 0
        total = integrand0
    elseif order == 1
        #total = integrand0 + (α/(4*π))*integrand1   
        total = integrand1
    elseif order == 2
        #total = integrand0 + (α/(4*π))*integrand1 + (α/(4*π))^2*integrand2
        total = integrand2
    end

    return total
end

#gluon---------------------------------------------------------------------------------------------

function Hg1_func(x)

    H0 = H0_func(x)
    H1 = H1_func(x)

    value = 1/4*CF*(
        (8*(-2*x + 2 + x^2)/x)*H0 - (4*(-2*x + 2 + x^2)/x)*H1
    )

    return value
end

function Hg2_func(x)

    logx = log(x)    
    log1mx = log(1-x)
    log1px = log(1+x)

    value = (
        (-32*(10 + 54*x + 27*x^2)*log1mx*reli(2, 1 - x))/(9*x) + (8*(2376 + 3192*x + 540*x^2 + 432*x^3 + 96*x^4 + 1296*log1mx + 1296*x*log1mx + 648*x^2*log1mx + 72*logx + 2448*x*logx + 72*x^2*logx)*
        reli(2, -x))/(81*x) + (8*(8424 - 3024*x - 468*x^2 - 432*x^3 + 720*log1mx - 720*x*log1mx + 360*x^2*log1mx - 4824*logx + 4752*x*logx - 2376*x^2*logx)*reli(2, x))/(81*x) + 
        (8*(-1296 - 1296*x - 648*x^2)*reli(3, (1 - x)/2))/(81*x) + (8*(1656 + 3240*x + 1620*x^2)*reli(3, 1 - x))/(81*x) + (8*(1224 - 3456*x + 576*x^2)*reli(3, -x))/(81*x) + 
        (8*(12600 - 4320*x + 8640*x^2)*reli(3, x))/(81*x) + (8*(1296 + 1296*x + 648*x^2)*reli(3, (2*x)/(-1 + x)))/(81*x) + (8*(-1296 - 1296*x - 648*x^2)*reli(3, x/(1 + x)))/(81*x) + 
        (8*(1296 + 1296*x + 648*x^2)*reli(3, (2*x)/(1 + x)))/(81*x) + (8*(-1296 - 1296*x - 648*x^2)*reli(3, (1 + x)/2))/(81*x) + 
        (8*(13497 - 162*π^2 - 7866*x - 264*π^2*x - 5601*x^2 + 216*π^2*x^2 - 1056*x^3 + 16*π^2*x^4 + 432*log2^3 + 432*x*log2^3 + 216*x^2*log2^3 - 108*π^2*log4 - 108*π^2*x*log4 - 54*π^2*x^2*log4 - 738*log1mx - 
            180*π^2*log1mx + 2214*x*log1mx - 24*π^2*x*log1mx + 18*x^2*log1mx - 96*π^2*x^2*log1mx - 648*x^3*log1mx - 648*log2^2*log1mx - 648*x*log2^2*log1mx - 324*x^2*log2^2*log1mx - 
            2538*log1mx^2 + 2304*x*log1mx^2 - 576*x^2*log1mx^2 + 216*x^3*log1mx^2 + 216*log8*log1mx^2 + 216*x*log8*log1mx^2 + 108*x^2*log8*log1mx^2 + 132*log1mx^3 - 564*x*log1mx^3 + 
            66*x^2*log1mx^3 + 2808*logx + 216*π^2*logx + 16434*x*logx - 768*π^2*x*logx + 2406*x^2*logx - 156*π^2*x^2*logx + 840*x^3*logx - 36*x^2*log1mx*logx + 1908*log1mx^2*logx - 
            1764*x*log1mx^2*logx + 558*x^2*log1mx^2*logx - 5724*logx^2 - 2736*x*logx^2 + 1575*x^2*logx^2 - 216*x^3*logx^2 - 48*x^4*logx^2 + 108*log1mx*logx^2 - 144*x*log1mx*logx^2 + 
            72*x^2*log1mx*logx^2 - 2160*logx^3 - 924*x*logx^3 - 1806*x^2*logx^3 + 612*π^2*log1px + 576*π^2*x*log1px + 288*π^2*x^2*log1px - 648*log2^2*log1px - 648*x*log2^2*log1px - 
            324*x^2*log2^2*log1px + 2376*logx*log1px + 3192*x*logx*log1px + 540*x^2*logx*log1px + 432*x^3*logx*log1px + 96*x^4*logx*log1px + 1332*logx^2*log1px + 1368*x*logx^2*log1px + 
            684*x^2*logx^2*log1px + 108*log(64)*log1px^2 + 108*x*log(64)*log1px^2 + 54*x^2*log(64)*log1px^2 - 72*(16 + 14*x + 7*x^2)*reli(3, 1 + x) - 9108*z3 + 3204*x*z3 - 6930*x^2*z3))/(81*x)
    )
    return value
end

function Hg_integrand(; x::Float64, α::Float64, order::Int64, fx::Float64, f1::Float64, nf::Int64)

    integrand0 = 0

    Hg1 = Hg1_func(x)
    integrand1 = (
        Hg1*fx
    )

    Hg2 = Hg2_func(x)
    integrand2 = (
        Hg2*fx
    )/4

    if order < 0
        total = 0
    elseif order == 0
        total = integrand0
    elseif order == 1
        #total = integrand0 + (α/(4*π))*integrand1
        total = integrand1
    elseif order == 2
        #total = integrand0 + (α/(4*π))*integrand1 + (α/(4*π))^2*integrand2
        total = integrand2
    end

    return total
end

f(x) = Hq_integrand(x=x, α=0.118, order=0, fx=x^2, f1=1.0, nf=5)
eps=10^(-12)
integral, error = quadgk(f, eps, 1-eps)
println(integral)