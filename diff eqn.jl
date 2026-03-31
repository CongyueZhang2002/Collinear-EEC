using DifferentialEquations
using ForwardDiff
using SpecialFunctions

function coupled_system!(dv, v, params, μ)
    Dq, Dg = v  # unpack state variables

    as = alpha_s_func(μ)/(4π)

    as2 = as*as
    as3 = as2*as

    #γqq0 = 5.55556
    #γqg0 = -2.33333
    #γgq0 = -1.55556
    #γgg0 = 11.7333

    #γqq1 = 32.4387
    #γqg1 = 19.7872
    #γgq1 = -44.8857
    #γgg1 = -79.9414

    #γqq2 = 50.2955
    #γqg2 = -130.6
    #γgq2 = -373.553
    #γgg2 = -675.797

    γqq0 = 8.088888864888888
    γqg0 = -1.5238095138095236
    γgq0 = -0.7111111084444444
    γgg0 = 18.847618990952377

    γqq1 = 45.46862898454472
    γqg1 = 7.98952645017585
    γgq1 = -23.929667885978443
    γgg1 = -25.427116808138297

    γqq2 = 99.29699567195216
    γqg2 = -166.85282742266804
    γgq2 = -393.74520800858716
    γgg2 = -905.5552656712619

    if nloops_αs == 1
        γqq = as*γqq0
        γqg = as*γqg0 
        γgq = as*γgq0 
        γgg = as*γgg0
    elseif nloops_αs == 2
        γqq = as*γqq0 + as2*γqq1
        γqg = as*γqg0 + as2*γqg1
        γgq = as*γgq0 + as2*γgq1
        γgg = as*γgg0 + as2*γgg1
    elseif nloops_αs == 3
        γqq = as*γqq0 + as2*γqq1 + as3*γqq2
        γqg = as*γqg0 + as2*γqg1 + as3*γqg2
        γgq = as*γgq0 + as2*γgq1 + as3*γgq2
        γgg = as*γgg0 + as2*γgg1 + as3*γgg2
    end

    dDq = -2/μ*(γqq*Dq+γgq*Dg)
    dDg = -2/μ*(γqg*Dq+γgg*Dg)

    dv[1] = dDq
    dv[2] = dDg
end

function test_gamma3(; kt, μ)

    Jq0 = 1.0;
    Jq1 = -16.444444444444443;
    Jq2 = -176.0183393572608;

    Jg0 = 1.0;
    Jg1 = -38.72;
    Jg2 = 181.0269818225895;

    as_0 = alpha_s_func(kt)/(4π)

    if nloops_αs == 1
        jet_q0 = Jq0
        jet_g0 = Jg0
    elseif nloops_αs == 2
        jet_q0 = Jq0 + as_0*Jq1
        jet_g0 = Jg0 + as_0*Jg1
    elseif nloops_αs == 3
        jet_q0 = Jq0 + as_0*Jq1 + as_0^2*Jq2
        jet_g0 = Jg0 + as_0*Jg1 + as_0^2*Jg2
    end

    dv0 = [jet_q0*FNP, jet_g0*FNP]

    μspan = (kt, μ)


    problem = ODEProblem(coupled_system!, dv0, μspan)
    solution = solve(problem)

    return solution.u[end][1], solution.u[end][2]
end

Jq1 = -16.444444444444443
Jq2 = -176.0183393572608
Jg1 = -38.72
Jg2 = 181.0269818225895

qL1 = -4.0
qL2 = -11.533333333333333
qL3 = -30.09925925925926
qL4 = -128.0129135802469
qL5 = -641.0409872153635
qL6 = -3511.5755904128437
qL7 = -20380.953209261206
qL8 = -123189.8134731188
qL9 = -767393.6244906444

qN1 = -16.444444444444443
qN2 = -112.70625665079928
qN3 = -728.9362049353601
qN4 = -4474.783338662404
qN5 = -29465.49297261155
qN6 = -200707.80484812864
qN7 = -1.3965377887472804*10^6
qN8 = -9.861081275849514*10^6
qN9 = -7.03851951087887*10^7

qNN2 = -176.0183393572608
qNN3 = -3494.5303089695185
qNN4 = -38488.004250181075
qNN5 = -376104.1656894854
qNN6 = -3.4269162265619566*10^6
qNN7 = -3.0074570254039057*10^7
qNN8 = -2.5761804466196218*10^8
qNN9 = -2.1696351203778615*10^9

function jet_analytic(; kt, μ, order)

    L = 2*log(μ/kt)
    as = alpha_s_func(μ)/(4π)

    if order == 0
        Jq = (
            1 + 
            as*(qL1*L) + 
            as^2*(L^2*qL2) + 
            as^3*(L^3*qL3) + 
            as^4*(L^4*qL4) + 
            as^5*(L^5*qL5) + 
            as^6*(L^6*qL6) + 
            as^7*(L^7*qL7) + 
            as^8*(L^8*qL8) + 
            as^9*(L^9*qL9)
        )
    elseif order == 1
        Jq = (
            1 + 
            as*(qL1*L + Jq1) + 
            as^2*(L^2*qL2 + L*qN2) + 
            as^3*(L^3*qL3 + L^2*qN3) + 
            as^4*(L^4*qL4 + L^3*qN4) + 
            as^5*(L^5*qL5 + L^4*qN5) + 
            as^6*(L^6*qL6 + L^5*qN6) + 
            as^7*(L^7*qL7 + L^6*qN7) + 
            as^8*(L^8*qL8 + L^7*qN8) + 
            as^9*(L^9*qL9 + L^8*qN9)
        )
    elseif order == 2
        Jq = (
            1 + 
            as*(qL1*L + Jq1) + 
            as^2*(L^2*qL2 + L*qN2 + Jq2) + 
            as^3*(L^3*qL3 + L^2*qN3 + L*qNN3) + 
            as^4*(L^4*qL4 + L^3*qN4 + L^2*qNN4) + 
            as^5*(L^5*qL5 + L^4*qN5 + L^3*qNN5) + 
            as^6*(L^6*qL6 + L^5*qN6 + L^4*qNN6) + 
            as^7*(L^7*qL7 + L^6*qN7 + L^5*qNN7) + 
            as^8*(L^8*qL8 + L^7*qN8 + L^6*qNN8) + 
            as^9*(L^9*qL9 + L^8*qN9 + L^7*qNN9)
        )
    end
    return Jq

end

function Diff_func(; z::Float64, Q::Float64, μf::Float64)

    Q2 = Q^2
    
    μ0 = 4.18

    N = 0.035
    a = 150
    b = 10.1

    Diffq0 = N*Q2*exp(-z*Q2/a)/(1+z*Q2/b)

    dv0 = [Diffq0,Diffq0]

    μspan = (μ0, μf)

    problem = ODEProblem(coupled_system!, dv0, μspan)
    solution = solve(problem)

    return solution.u[end][1], solution.u[end][2]
end

function jet_b_analytic(; b, μ, order, a1, a2)

    global nloops_αs = order + 1
    bstar = b/sqrt(1 + (b/b0)^2)

    init_splitting_globals(b=bstar, a1=a1, a2=a2)

    Lb = 2*log(bstar*μ/b0)
    as = alpha_s_func(μ)/(4π)

    return Jq_NP_func(a_s=as, Lb=Lb, order=order)
end

function jet_b_numerical(; b, μ, order, a1, a2)

    global nloops_αs = order + 1
    bstar = b/sqrt(1 + (b/b0)^2)
    μ0 = b0/bstar

    init_splitting_globals(b=bstar, a1=a1, a2=a2)

    return test_gamma3(kt = μ0, μ=μ)[1]
end

function jet_kt_analytic(; kt, μ, order, a1, a2)
    integrand(b) = b*besselj0(b*kt)/(2*pi) * jet_b_analytic(b=b, μ=μ, order=order, a1=a1, a2=a2)
    return quadgk(integrand,0.001,5.0,rtol=0.005)[1]
end

function jet_kt_numerical(; kt, μ, order, a1, a2)
    integrand(b) = b*besselj0(b*kt)/(2*pi) * jet_b_numerical(b=b, μ=μ, order=order, a1=a1, a2=a2)
    return quadgk(integrand,0.001,5.0,rtol=0.005)[1]
end