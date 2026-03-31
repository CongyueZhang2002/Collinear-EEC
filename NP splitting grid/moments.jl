using PolyLog
using QuadGK
include("../core/constants.jl")

function splittings_func(; x::Float64, type::String)

    x2 = x^2
    x3 = x^3
    x4 = x^4

    logx = log(x)
    log1mx = log(1 - x)
    log1px = log(1 + x)

    #LO

    if type == "Pqq0"
        value = (
            ((-8*(1 + x))/3), 
            (4), 
            (16/3)
        )

    elseif type == "Pqg0"
        value = (5*(2 - 4*x + 4*x2))

    elseif type == "Pgq0"
        value = ((4*(-4 + 4/x + 2*x))/3)

    elseif type == "Pgg0"
        value = (
            (12*(-2 + x^(-1) + x - x2)), 
            (23/3), 
            (12)
        )

    #NLO

    elseif type == "Pqq1"
        value = (
            ((8*(-200 + 6*logx^2*x*(5 + 4*x)*(4 + 5*x) + x*(-903 + 2*π^2*(5 + x*(9 + 5*x)) + x*(-662 + x*(601 + 560*x))) + 3*logx*x*(-((1 + x)*(197 + 16*log1mx*(1 + x) + x*(301 + 80*x))) + 4*log1px*(1 + x2)) + 
                12*(x + x3)*reli(2,-x)))/(27*x*(1 + x))),
            ((2*(159 + 44*π^2 - 72*z3))/27),
            ((-16*(-151 + 21*(-5 + logx)*logx + 9*π^2))/27),
            ((256*logx)/9)
        )

    elseif type == "Pqg1"
        value = ((-20*(60 + 285*x2 + x*(83 - 178*x - 3*logx^2*(7 + 58*x - 8*x2) + 15*log1mx^2*(-1 + 2*x - 2*x2) + 160*x2 + 2*π^2*(2 + 5*x + 4*x2) + log1mx*(-33 + 6*x - 6*x2 + 36*logx*(1 - 2*x + 2*x2)) + 
                6*logx*(16 + 37*x + 12*x2 + 9*log1px*(1 + 2*x + 2*x2))) - 534*x3 + 6*x*(9*(1 + 2*x + 2*x2)*reli(2,-x) + 20*(1 - 2*x + 2*x2)*reli(2,x))))/(9*x))

    elseif type == "Pgq1"
        value = ((16*(17 + x*(43 + 6*π^2 + 18*x) + 5*log1mx^2*(-2 + 2*x - x2) - 9*x2 + log1mx*(-10*x^2 + 38*logx*(2 - 2*x + x2)) - logx^2*(2*x*(7 + x) + 9*(4 + 3*x2)) - 44*x3 + 
                logx*(-54 + 2*x*(20 + x) + 81*x2 + 18*log1px*(2 + 2*x + x2) + 24*x3) + 18*(2 + 2*x + x2)*reli(2,-x) + 80*reli(2,x) + 40*(-2*x + x2)*reli(2,x)))/(9*x))

    elseif type == "Pgg1"
        value = (
            ((1970*x - 2*(230 + 1278*x2 + 3703*x3 + 1210*x4) - 12*(3*logx^2*(61*x + 2*(54 + 88*x2 + 71*x3 - 27*x4)) + 
                logx*(614 + 935*x + 252*x2 - 324*log1px*(1 + x + x2)^2 + 185*x3 + 254*x4 + 324*log1mx*(-1 + x + x2 + x4)) - 27*π^2*(3*x + 2*(2*x2 + x3 + x4))) + 3888*(1 + x + x2)^2*reli(2,-x))/(27*x*(1 + x))),
            ((4*(32 + 81*z3))/3),
            (604/3 + 4*(46 - 27*logx)*logx - 12*π^2),
            (144*logx)
        )

        #NNLO

    elseif type == "Pqq2"
        value = (
            (log1mx*(-19.595999999999886 + logx*(117.58059999999995 + logx*(-376.61499999999995 - 115.50999999999999*x) - 558.7775000000001*x) - 372.025*x) + 
                log1mx^2*(-4.305250000000002 + 4.305250000000002*x) + log1mx^3*(-29.630000000000003 + 29.630000000000003*x) + 
                (2359.838999999997 + logx*(1620.344999999998 + x*(3696.609999999995 + (-4372.124999999996 - 657.3474999999993*x)*x) + 
                logx*(-71.11099999999999 + (2240.5947499999997 - 2369.71975*x)*x + 
                logx*(-142.22199999999964 + x*(304.0728999999991 + logx*(47.41829 - 45.838049999999996*x) + x*(-385.722249999999 + 196.61249999999953*x))))) + 
                x*(108.24250000000075 + x*(-9385.124999999989 + x*(6315.922499999994 + x*(-29.295000000000755 + x*(-534.1049999999996 + 199.83499999999984*x))))))/x),
            (454.2215),
            (239.20192500000007 + logx*(98.76525000000001 + 29.629499999999997*logx)*x)
        )

    elseif type == "Pqg2"
        value = ((5636.035000000001 + logx^4*(211.64000000000004 - 4427.5*x)*x + logx^3*(-320.00000000000006 + (2345.625 - 8345.0*x)*x) + 
                logx^2*(-349.6296296296297 + x*(14117.166666666668 - 333.33333333333337*x + 333.33333333333337*x2)) + 
                logx*(3299.770000000001 + x*(25369.305555555555 - 444.4444444444445*x + 444.4444444444445*x2 + 
                log1mx*(15460.333333333338 - 39.825*log1mx - 18932.11666666667*x + 666.6666666666667*x2))) + 
                x*(610.5 + 3689.1444444444455*x + 2712.2055555555567*x2 + log1mx*(1434.530555555556 - 444.4444444444445*x + 444.4444444444445*x2 + 
                    log1mx*(595.8666666666669 + log1mx*(101.85185185185185 + 18.51851851851852*log1mx) - 333.33333333333337*x + 333.33333333333337*x2)) - 
                    11426.499999999998*x3 - 506.32500000000005*x4))/x)

    elseif type == "Pgq2"
        value = ((5833.720499999999 + logx^2*(1541.1492592592592 - 1716.2125200000003*x) + logx^4*(256.0 - 30.061*x) + 
                logx*(5486.15 + (-2455.5864 + (-1025.3599999999997 - 71.23*log1mx)*log1mx)*x) + logx^3*(1316.3456790123455 + x*(-369.38000000000005 + 522.0675*x)) + 
                x*(-4083.57 + log1mx*(-237.17000000000002 + log1mx*(-155.93246913580248 + log1mx*(24.197530864197528 + 4.938271604938271*log1mx))) - 3238.049999999999*x + 
                    4487.05*x2 + 363.24*x3 + 128.19*x4))/x)

    elseif type == "Pgg2"
        value = (
            (438.8550000000001*log1mx^2*logx + log1mx*(-1990.25 + logx*(3270.8749999999995 - 23049.245*logx + 2413.0474999999997*x)) + 
                (10242.3095 + logx*(10319.23 + 10540.167*x + logx*(4996.2845 + 14822.7*x + logx*(3416.8885 + logx*(576.0 + 282.415*x) + x*(3931.275 + 8111.7*x)))) + 
                x*(-29174.002500000002 + x*(12603.75 + x*(29134.225 + x*(-42051.125 + 14921.654999999999*x)))))/x),
            (1943.4259999999997),
            (538.2157499999998)
                )
    end

    return value
end

function integrand_func(; y::Float64, f::Function, type::String)

    Pqq0Delta_at_1 = 4
    Pqq0D0_at_1 = 16/3

    Pgg0Delta_at_1 = 23/3
    Pgg0D0_at_1 = 12

    Pqq1Delta_at_1 = 37.534407157069694
    Pqq1D0_at_1 = 36.843591342338236
    Pqq1D1_at_1 = 0

    Pgg1Delta_at_1 = 172.48881220790284
    Pgg1D0_at_1 = 82.89808052026105
    Pgg1D1_at_1 = 0

    Pqq2Delta_at_1 = 454.2215
    Pqq2D0_at_1 = 239.201925
 
    Pgg2Delta_at_1 = 1943.426
    Pgg2D0_at_1 = 538.21575

    factor = f(y)
    factor_at_1 = f(1)

    #LO

    if type == "Pqq0" 

        Pqq0Reg, Pqq0Delta, Pqq0D0 = splittings_func(x=y, type="Pqq0")
        value = (
            Pqq0Reg * factor
            + Pqq0Delta_at_1 * factor_at_1
            + (Pqq0D0 * factor - Pqq0D0_at_1 * factor_at_1) * 1/(1-y) 
        )

    elseif type == "Pqg0"

        Pqg0 = splittings_func(x=y, type="Pqg0")    
        value = (
            Pqg0 * factor
        )

    elseif type == "Pgq0"

        Pgq0 = splittings_func(x=y, type="Pgq0")        
        value = (
            Pgq0 * factor
        )

    elseif type == "Pgg0"

        Pgg0Reg, Pgg0Delta, Pgg0D0 = splittings_func(x=y, type="Pgg0")
        value = (
            Pgg0Reg * factor
            + Pgg0Delta_at_1 * factor_at_1
            + (Pgg0D0 * factor - Pgg0D0_at_1 * factor_at_1) * 1/(1-y)
        )

    # NLO

    elseif type == "Pqq1"

        Pqq1Reg, Pqq1Delta, Pqq1D0, Pqq1D1 = splittings_func(x=y, type="Pqq1")
        value = (
            Pqq1Reg * factor
            + Pqq1Delta_at_1 * factor_at_1
            + (Pqq1D0 * factor - Pqq1D0_at_1 * factor_at_1) * 1/(1-y)
            + (Pqq1D1 * factor - Pqq1D1_at_1 * factor_at_1) * log(1-y)/(1-y)
        )

    elseif type == "Pqg1"

        Pqg1 = splittings_func(x=y, type="Pqg1")
        value = (
            Pqg1 * factor
        )

    elseif type == "Pgq1"

        Pgq1 = splittings_func(x=y, type="Pgq1")
        value = (
            Pgq1 * factor
        )

    elseif type == "Pgg1"

        Pgg1Reg, Pgg1Delta, Pgg1D0, Pgg1D1 = splittings_func(x=y, type="Pgg1")        
        value = (
            Pgg1Reg * factor
            + Pgg1Delta_at_1 * factor_at_1
            + (Pgg1D0 * factor - Pgg1D0_at_1 * factor_at_1) * 1/(1-y)
            + (Pgg1D1 * factor - Pgg1D1_at_1 * factor_at_1) * log(1-y)/(1-y)
        )

    #NNLO

    elseif type == "Pqq2"

        Pqq2Reg, Pqq2Delta, Pqq2D0 = splittings_func(x=y, type="Pqq2")
        value = (
            Pqq2Reg * factor
            + Pqq2Delta_at_1 * factor_at_1
            + (Pqq2D0 * factor - Pqq2D0_at_1 * factor_at_1) * 1/(1-y)
        )

    elseif type == "Pqg2"

        Pqg2 = splittings_func(x=y, type="Pqg2")
        value = (
            Pqg2 * factor
        )

    elseif type == "Pgq2"

        Pgq2 = splittings_func(x=y, type="Pgq2")
        value = (
            Pgq2 * factor
        )

    elseif type == "Pgg2"

        Pgg2Reg, Pgg2Delta, Pgg2D0 = splittings_func(x=y, type="Pgg2")
        value = (
            Pgg2Reg * factor
            + Pgg2Delta_at_1 * factor_at_1
            + (Pgg2D0 * factor - Pgg2D0_at_1 * factor_at_1) * 1/(1-y)
        )

    end

    return value
end

function moments_func(fNP::Function)

    factor(y) = -fNP(y)*y^2

    γqq0_integrand(y) = integrand_func(y=y, f=factor, type="Pqq0")
    γqg0_integrand(y) = integrand_func(y=y, f=factor, type="Pqg0")
    γgq0_integrand(y) = integrand_func(y=y, f=factor, type="Pgq0")
    γgg0_integrand(y) = integrand_func(y=y, f=factor, type="Pgg0")

    γqq1_integrand(y) = integrand_func(y=y, f=factor, type="Pqq1")
    γqg1_integrand(y) = integrand_func(y=y, f=factor, type="Pqg1")
    γgq1_integrand(y) = integrand_func(y=y, f=factor, type="Pgq1")
    γgg1_integrand(y) = integrand_func(y=y, f=factor, type="Pgg1")

    γqq2_integrand(y) = integrand_func(y=y, f=factor, type="Pqq2")
    γqg2_integrand(y) = integrand_func(y=y, f=factor, type="Pqg2")
    γgq2_integrand(y) = integrand_func(y=y, f=factor, type="Pgq2")
    γgg2_integrand(y) = integrand_func(y=y, f=factor, type="Pgg2")

    factorln(y) = -fNP(y)*y^2*log(y)

    dγqq0_integrand(y) = integrand_func(y=y, f=factorln, type="Pqq0")
    dγqg0_integrand(y) = integrand_func(y=y, f=factorln, type="Pqg0")
    dγgq0_integrand(y) = integrand_func(y=y, f=factorln, type="Pgq0")
    dγgg0_integrand(y) = integrand_func(y=y, f=factorln, type="Pgg0")

    dγqq1_integrand(y) = integrand_func(y=y, f=factorln, type="Pqq1")
    dγqg1_integrand(y) = integrand_func(y=y, f=factorln, type="Pqg1")
    dγgq1_integrand(y) = integrand_func(y=y, f=factorln, type="Pgq1")
    dγgg1_integrand(y) = integrand_func(y=y, f=factorln, type="Pgg1")

    factorln2(y) = -fNP(y)*y^2*log(y)^2
    
    ddγqq0_integrand(y) = integrand_func(y=y, f=factorln2, type="Pqq0")
    ddγqg0_integrand(y) = integrand_func(y=y, f=factorln2, type="Pqg0")
    ddγgq0_integrand(y) = integrand_func(y=y, f=factorln2, type="Pgq0")
    ddγgg0_integrand(y) = integrand_func(y=y, f=factorln2, type="Pgg0")

    eps = 1e-8  

    γqq0 = quadgk(γqq0_integrand, eps, 1-eps)[1]
    γqg0 = quadgk(γqg0_integrand, eps, 1-eps)[1]
    γgq0 = quadgk(γgq0_integrand, eps, 1-eps)[1]
    γgg0 = quadgk(γgg0_integrand, eps, 1-eps)[1]

    dγqq0 = quadgk(dγqq0_integrand, eps, 1-eps)[1]
    dγqg0 = quadgk(dγqg0_integrand, eps, 1-eps)[1]
    dγgq0 = quadgk(dγgq0_integrand, eps, 1-eps)[1]
    dγgg0 = quadgk(dγgg0_integrand, eps, 1-eps)[1]

    ddγqq0 = quadgk(ddγqq0_integrand, eps, 1-eps)[1]
    ddγqg0 = quadgk(ddγqg0_integrand, eps, 1-eps)[1]
    ddγgq0 = quadgk(ddγgq0_integrand, eps, 1-eps)[1]
    ddγgg0 = quadgk(ddγgg0_integrand, eps, 1-eps)[1]

    γqq1 = quadgk(γqq1_integrand, eps, 1-eps)[1]
    γqg1 = quadgk(γqg1_integrand, eps, 1-eps)[1]
    γgq1 = quadgk(γgq1_integrand, eps, 1-eps)[1]
    γgg1 = quadgk(γgg1_integrand, eps, 1-eps)[1]

    dγqq1 = quadgk(dγqq1_integrand, eps, 1-eps)[1]
    dγqg1 = quadgk(dγqg1_integrand, eps, 1-eps)[1]
    dγgq1 = quadgk(dγgq1_integrand, eps, 1-eps)[1]
    dγgg1 = quadgk(dγgg1_integrand, eps, 1-eps)[1]

    γqq2 = quadgk(γqq2_integrand, eps, 1-eps)[1]
    γqg2 = quadgk(γqg2_integrand, eps, 1-eps)[1]
    γgq2 = quadgk(γgq2_integrand, eps, 1-eps)[1]
    γgg2 = quadgk(γgg2_integrand, eps, 1-eps)[1]

    return (
        γqq0, γqg0, γgq0, γgg0,
        dγqq0, dγqg0, dγgq0, dγgg0,
        ddγqq0, ddγqg0, ddγgq0, ddγgg0,
        γqq1, γqg1, γgq1, γgg1,
        dγqq1, dγqg1, dγgq1, dγgg1,
        γqq2, γqg2, γgq2, γgg2
    )

end

function buildgrid_func(; a1b::Float64, a2::Float64)

    fNP(y) = exp(-(a1b/y)^a2)
    
    return moments_func(fNP)
end

#a1 = 2.5
#a2 = 1.1
#b = 5.0
#fNP(y) = exp(-(a1*b/y)^a2)
#fNP(y) = 1
#display(moments_func(fNP))