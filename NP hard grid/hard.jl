#https://www.sciencedirect.com/science/article/pii/S0550321306004330
#timelike singlet splitting functions

using QuadGK
using PolyLog
include("harmonic polylogs.jl")
#include("constants.jl")

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

    logx2 = logx^2  
    log1mx2 = log1mx^2
    log1px2 = log1px^2

    x2 = x^2
    x3 = x^3
    x4 = x^4
    x5 = x^5
    x6 = x^6

    value = (
        1.7777777777777777*(41.11111111111111 - 20*logx + 20*log1px*logx - 40*logx2 + 251.51866124483652*x - 9.389151387171156*log1px*x - 
        0.6931471805599452*log1px2*x + 188.56849445704398*logx*x + 59*log1px*logx*x - 9.875*logx2*x - 2*log1px*logx2*x - 23.833333333333332*logx^3*x + 
        16.24638151918887*x2 + 9.389151387171156*log1px*x2 + 0.6931471805599452*log1px2*x2 + 119.85836259074044*logx*x2 + 37*log1px*logx*x2 + 48.625*logx2*x2 + 
        2*log1px*logx2*x2 - 25.333333333333332*logx^3*x2 - 31.613594924924225*x3 - 19.258755788260515*log1px*x3 - 0.6931471805599452*log1px2*x3 - 
        104.38713911826387*logx*x3 - 39*log1px*logx*x3 - 27.375*logx2*x3 -log1px*logx2*x3 + 38.5*logx^3*x3 + 
        log1mx^3*(4.333333333333333*x + 3.6666666666666665*x2 - 3.6666666666666665*x3 - 4.333333333333333*x4) + 199.75188634579806*x4 + 19.258755788260515*log1px*x4 + 
        0.6931471805599452*log1px2*x4 - 50.43060431790076*logx*x4 - 59*log1px*logx*x4 - 34.875*logx2*x4 +log1px*logx2*x4 + 37*logx^3*x4 + 
        log1mx2*(-10+ (-25.25 - 0.3333333333333333*log8 - 18*logx)*x + (7.75 + 0.3333333333333333*log8 - 16*logx)*x2 + 35.25*x3 - 0.3333333333333333*log8*x3 + 
            16*logx*x3 + 2.25*x4 + 0.3333333333333333*log8*x4 + 18*logx*x4 - 10*x5) + 39.11111111111111*x5 - 45.33333333333333*logx*x5 - 20*log1px*logx*x5 + 
        10*logx2*x5 + log1mx*(-10+ (113.80284494525438 - 20*logx + 1.5*logx2)*x + (-54.947929216278474 + 52*logx + 1.5*logx2)*x2 - 192.45075212068605*x3 + 
            44*logx*x3 + 0.5*logx2*x3 + 27.798473717781093*x4 - 28*logx*x4 + 0.5*logx2*x4 + 50*x5) + 3.2898681336964524*x6 + 2*log1px*logx*x6 -logx2*x6 + 
        49*log1mx*x*(1+ x)*(-0.3877551020408163 +x2)*reli(2, 1-x) - 2*(-10+ (-29.5 +log1mx + 2*logx)*x + (-18.5 -log1mx)*x2 + 19.5*x3 + 
            log1mx*x3 + 3*logx*x3 + 29.5*x4 -log1mx*x4 +logx*x4 + 10*x5 -x6)*reli(2, -1*x) + 
        10*(4+ (-3.65 +log1mx + 4.8*logx)*x + (-8.25 +log1mx + 4.8*logx)*x2 + 0.55*x3 +log1mx*x3 + 4.6*logx*x3 + 3.15*x4 +log1mx*x4 + 4.6*logx*x4 + 
            2*x5)*reli(2, x) + 17*(1*x + 1.2352941176470589*x2 - 3*x3 - 2.764705882352941*x4)*reli(3, 1-x) + 
        2*(1*x -x2 +x3 -x4)*reli(3, 0.5 - 0.5*x) + 2*(1*x -x2 +x3 -x4)*reli(3, 0.5 + 0.5*x) + 
        6*(1*x + 0.3333333333333333*x2 + 2*x3 + 0.6666666666666666*x4)*reli(3, -1*x) + 4*(1*x + 1.5*x2 - 49.5*x3 - 49*x4)*reli(3, x) - 
        2*(1*x -x2 +x3 -x4)*reli(3, (2*x)/(-1+ x)) + 2*(1*x -x2 +x3 -x4)*reli(3, x/(1+ x)) - 
        2*(1*x -x2 +x3 -x4)*reli(3, (2*x)/(1+ x)) + x*(2- 2*x + 6*x2 - 6*x3)*reli(3, 1+ x)))/(x*(-1+ x2)
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

function Hq_integrand(; x::Float64, α::Float64, order::Int64, fNP::Function, f1::Float64, nf::Int64)

    fx = fNP(x)

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

    logx2 = logx^2
    log1mx2 = log1mx^2
    log1px2 = log1px^2
    
    x2 = x^2
    x3 = x^3
    x4 = x^4
    x5 = x^5
    x6 = x^6

    value = (
        -37.92776745560282 - 279.0975155767979*log1mx - 250.66666666666666*log1mx2 + 13.037037037037036*log1mx^3 + 565.8137620195251*log1px + 
        44.36141955583649*log1px2 + 21.333333333333332*log1mx2*log8 + 487.8848938899063*logx + 188.44444444444443*log1mx2*logx + 234.66666666666666*log1px*logx - 
        565.3333333333333*logx2 + 10.666666666666666*log1mx*logx2 + 131.55555555555554*log1px*logx2 - 213.33333333333331*logx^3 - 785.5796249578439*x + 
        164.52305593628256*log1mx*x + 227.55555555555554*log1mx2*x - 55.7037037037037*log1mx^3*x + 530.7218352600963*log1px*x + 44.36141955583649*log1px2*x + 
        21.333333333333332*log1mx2*log8*x + 874.4833402432962*logx*x - 174.2222222222222*log1mx2*logx*x + 315.25925925925924*log1px*logx*x - 
        270.22222222222223*logx2*x - 14.222222222222221*log1mx*logx2*x + 135.11111111111111*log1px*logx2*x - 91.25925925925925*logx^3*x - 1231.2418782743352*x2 - 
        107.17519002608154*log1mx*x2 - 56.888888888888886*log1mx2*x2 + 6.518518518518518*log1mx^3*x2 + 265.36091763004816*log1px*x2 + 22.180709777918246*log1px2*x2 + 
        10.666666666666666*log1mx2*log8*x2 + 85.5646136721047*logx*x2 - 3.5555555555555554*log1mx*logx*x2 + 55.11111111111111*log1mx2*logx*x2 + 
        53.33333333333333*log1px*logx*x2 + 155.55555555555554*logx2*x2 + 7.111111111111111*log1mx*logx2*x2 + 67.55555555555556*log1px*logx2*x2 - 
        178.37037037037035*logx^3*x2 - 104.29629629629629*x3 - 64*log1mx*x3 + 21.333333333333332*log1mx2*x3 + 82.96296296296296*logx*x3 + 
        42.666666666666664*log1px*logx*x3 - 21.333333333333332*logx2*x3 + 15.596411893079479*x4 + 9.481481481481481*log1px*logx*x4 - 4.7407407407407405*logx2*x4 - 
        192*log1mx*(0.18518518518518517 +x + 0.5*x2)*reli(2, 1-x) + 
        128*(1.8333333333333333 + 2.462962962962963*x + logx*(0.05555555555555555 + 1.8888888888888888*x + 0.05555555555555555*x2) + log1mx*(1+x + 0.5*x2) + 
            0.4166666666666667*x2 + 0.3333333333333333*x3 + 0.07407407407407407*x4)*reli(2, -1*x) - 
        71.11111111111111*(-11.7 + 4.2*x + log1mx*(-1+x - 0.5*x2) + 0.65*x2 + logx*(6.7 - 6.6*x + 3.3*x2) + 0.6*x3)*reli(2, x) + 
        320*(0.5111111111111111 +x + 0.5*x2)*reli(3, 1-x) - 128*(1+x + 0.5*x2)*reli(3, 0.5 - 0.5*x) - 
        128*(1+x + 0.5*x2)*reli(3, 0.5 + 0.5*x) - 341.3333333333333*(-0.3541666666666667 +x - 0.16666666666666666*x2)*reli(3, -1*x) - 
        426.66666666666663*(-2.9166666666666665 +x - 2*x2)*reli(3, x) + 128*(1+x + 0.5*x2)*reli(3, (2*x)/(-1+ x)) - 
        128*(1+x + 0.5*x2)*reli(3, x/(1+ x)) + 128*(1+x + 0.5*x2)*reli(3, (2*x)/(1+ x)) - 
        99.55555555555554*(1.1428571428571428 +x + 0.5*x2)*reli(3, 1+ x)
        )/x
    return value
end

function Hg_integrand(; x::Float64, α::Float64, order::Int64, fNP::Function, f1::Float64, nf::Int64)

    fx = fNP(x)

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

function hard_moments_func(fNP::Function)

    factor(x) = fNP(x)*x^2
    f1 = factor(1.0)

    hq0_integrand(x) = Hq_integrand(x=x, α=0.118, order=0, fNP=factor, f1=f1, nf=5)
    hq1_integrand(x) = Hq_integrand(x=x, α=0.118, order=1, fNP=factor, f1=f1, nf=5)
    hq2_integrand(x) = Hq_integrand(x=x, α=0.118, order=2, fNP=factor, f1=f1, nf=5)

    hg0_integrand(x) = Hg_integrand(x=x, α=0.118, order=0, fNP=factor, f1=f1, nf=5)
    hg1_integrand(x) = Hg_integrand(x=x, α=0.118, order=1, fNP=factor, f1=f1, nf=5)
    hg2_integrand(x) = Hg_integrand(x=x, α=0.118, order=2, fNP=factor, f1=f1, nf=5)

    eps=1e-9

    hq0 = quadgk(hq0_integrand, eps, 1-eps)[1]
    hq1 = quadgk(hq1_integrand, eps, 1-eps)[1]
    hq2 = quadgk(hq2_integrand, eps, 1-eps)[1]

    hg0 = quadgk(hg0_integrand, eps, 1-eps)[1]
    hg1 = quadgk(hg1_integrand, eps, 1-eps)[1]
    hg2 = quadgk(hg2_integrand, eps, 1-eps)[1]

    return (
        hq0, hq1, hq2,
        hg0, hg1, hg2
    )

end

function hard_buildgrid_func(; a1b::Float64, a2::Float64)

    fNP(x) = exp(-(a1b/x)^a2)
    
    return hard_moments_func(fNP)
end