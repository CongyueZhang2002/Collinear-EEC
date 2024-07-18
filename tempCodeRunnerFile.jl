using Distributed
using HCubature
using SpecialFunctions

@everywhere begin
    
    using HCubature
    using SpecialFunctions

    include("strong coupling\\constants.jl")
    include("strong coupling\\alpha_s.jl")
    include("anomalous dim\\cusp\\cusp.jl")
    include("anomalous dim\\non-cusp\\gammaH.jl")
    include("anomalous dim\\non-cusp\\gammaB.jl")
    include("anomalous dim\\non-cusp\\gammaS.jl")
    include("anomalous dim\\non-cusp\\gammav.jl")
    include("anomalous dim\\non-cusp\\Integration.jl")    
    include("hard\\hard.jl")
    include("jet\\jet.jl")
    include("soft\\soft.jl")
    include("NP.jl")
    include("renormalon.jl")

# Resummation Accuracy, order = # of loops
    function ResAcc(XLL)

        if     XLL == "NLL"
            FO_order = 0
            RES_order = 1
        elseif XLL == "NLL'"
            FO_order = 1
            RES_order = 1
        elseif XLL == "NNLL"
            FO_order = 1
            RES_order = 2
        elseif XLL == "NNLL'"
            FO_order = 2
            RES_order = 2
        elseif XLL == "N3LL"
            FO_order = 2
            RES_order = 3            
        elseif XLL == "N3LL'"
            FO_order = 3
            RES_order = 3
        elseif XLL == "N4LL"
            FO_order = 3
            RES_order = 4
        end

        return FO_order, RES_order
    end

# Analytic

    function Integral_H_analytic(; Q::Float64, μf::Float64, μH::Float64, 
                                   αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

        γH0 = γH0_tilde(nf) 
        γH1 = γH1_tilde(nf)
        γH2 = γH2_tilde(nf)
        γH3 = γH3_tilde(nf)

        total = (
            - 4*KΓ_integral(; μi=μH, μf=μf, scale=Q, αs_ini=αs_ini, μ_ini=μ_ini, order=order+1, nf=nf)
            + η_integral(; μi=μH, μf=μf, 
                        F0=γH0, F1=γH1, F2=γH2, F3=γH3, F4=0.0,
                        αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
        )
    
        return total
    end

    function Integral_J_analytic(; νJ::Float64, μJ::Float64, μf::Float64, Q::Float64,
                                   αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

        γJ0 = γB0_tilde(nf) 
        γJ1 = γB1_tilde(nf) 
        γJ2 = γB2_tilde(nf) 
        γJ3 = γB3_tilde(nf) 

        Γ0 = Γ0_func(nf)
        Γ1 = Γ1_func(nf)   
        Γ2 = Γ2_func(nf)   
        Γ3 = Γ3_func(nf)
        Γ4 = Γ4_func(nf)

        total = (
          η_integral(; μi=μJ, μf=μf, 
                F0=Γ0, F1=Γ1, F2=Γ2, F3=Γ3, F4=Γ4,
                αs_ini=αs_ini, μ_ini=μ_ini, order=order+1, nf=nf) * 2*log(νJ/Q)
        + η_integral(; μi=μJ, μf=μf, 
                F0=γJ0, F1=γJ1, F2=γJ2, F3=γJ3, F4=0.0,
                αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
        )

        return total
    end

    function Integral_ν_analytic(; b::Float64, μf::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64, bmax::Float64)

        bstar = b/(1+(b/bmax)^2)^0.5
    
        μ0 = b0/bstar
    
        αs_μ0 = alpha_s_func(μf=μ0, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
    
        γν_FO = γν_func(b=b, μ0=μ0, αs=αs_μ0, order=order, nf=nf, bmax=bmax)
        
        Γ0 = Γ0_func(nf)
        Γ1 = Γ1_func(nf)   
        Γ2 = Γ2_func(nf)   
        Γ3 = Γ3_func(nf)
        Γ4 = Γ4_func(nf)   
     
        total = (
            γν_FO
            -4 * η_integral(; μi=μ0, μf=μf, 
                           F0=Γ0, F1=Γ1, F2=Γ2, F3=Γ3, F4=Γ4,
                           αs_ini=αs_ini, μ_ini=μ_ini, order=order+1, nf=nf)
        )
    
        return total
    
    end

    function Integrand(;b::Float64, Q::Float64, z::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String,
        μJ_ratio::Float64, νJ_ratio::Float64, μH_ratio::Float64,
        μS_ratio::Float64, νS_ratio::Float64, bmax_ratio::Float64,
        parameters::Vector{Float64})

        FO_order, RES_order = ResAcc(XLL)

        if b < 1.0e-9
            return 0.0
        end

        J0 = besselj0(b*Q*(1-z)^0.5)

        bmax = bmax_ratio * b0
        bstar = b/(1+(b/bmax)^2)^0.5

        μJ = μJ_ratio * b0/bstar
        νJ = νJ_ratio * Q
        μS = μS_ratio * b0/bstar
        νS = νS_ratio * b0
        μH = μH_ratio * Q
        
        μf = μS # μf is the μ

        αs_J = alpha_s_func(μf=μJ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
        αs_S = alpha_s_func(μf=μS, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
 
        J = J_func(b=b, μJ=μJ, νdQ=νJ/Q, αs=αs_J, order=FO_order, nf=nf)
        S = S_func(b=b, μS=μS, νS=νS, αs=αs_S, order=FO_order, nf=nf)

        γH_Integral = Integral_H_analytic(Q=Q, μf=μf, μH=μH, 
                                αs_ini=αs_ini, μ_ini=μ_ini, order=RES_order, nf=nf)

        γJ_Integral = Integral_J_analytic(νJ=νJ, μJ=μJ, μf=μf, Q=Q,
                                αs_ini=αs_ini, μ_ini=μ_ini, order=RES_order, nf=nf)
        
        #γS_Integral = Integral_S(νS=νS, μf=μf, μS=μS, 
        #                        αs_ini=αs_ini, μ_ini=μ_ini, order=RES_order, nf=nf)

        γν_Integral = Integral_ν_analytic(b=b, μf=μf, 
                                αs_ini=αs_ini, μ_ini=μ_ini, order=RES_order, nf=nf, bmax=bmax)

        NP_Sudakov = NP(b=b, parameters=parameters)

        total = Q^2/4*b*J0*J*J*S*exp(γH_Integral+2*γJ_Integral)*(νJ/νS)^γν_Integral*NP_Sudakov

        return total
    end

    function Hard_Part(; Q::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, μH_ratio::Float64, XLL::String)
        
        FO_order, RES_order = ResAcc(XLL)

        μH = μH_ratio * Q
        αs_H = alpha_s_func(μf=μH, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
        H = H_func(μH=μH, Q=Q, αs=αs_H, order=FO_order, nf=nf)

        return H
    end

    function R_ratio(; Q::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64)
        αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4, nf=nf)
        r = 1 + αs/π + (αs/(4π))^2*(CF*CA*(123/2-44*z3) + CF*TF*nf*(-22+16*z3) - CF^2*3/2)
        return r
    end
    print(R_ratio(Q=29.0, μ_ini=91.2, αs_ini=0.118, nf=5))
    function sigma_z(; Q::Float64, z::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String,
                    μH_ratio::Float64, μJ_ratio::Float64, νJ_ratio::Float64, 
                    μS_ratio::Float64, νS_ratio::Float64, bmax_ratio::Float64,
                    parameters::Vector{Float64}, 
                    Ω1::Float64, μ_ren_ratio::Float64)

        f(x) = Integrand(b=x[1], Q=Q, z=z, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL,
                        μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, μS_ratio=μS_ratio, 
                        νS_ratio=νS_ratio, μH_ratio=μH_ratio, bmax_ratio=bmax_ratio,
                        parameters=parameters)               

        integral, error = hcubature(f, [0.0], [10.0])
        
        hard = Hard_Part(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, μH_ratio=μH_ratio, XLL=XLL)

        μ_ren = μ_ren_ratio * Q
        renormalon = renormalon_MS_func(z=z, Q=Q, μ_ren=μ_ren, μ_ini=μ_ini, αs_ini=αs_ini, order=2, nf=nf, Ω1=Ω1) 

        r = R_ratio(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf)

        total = (hard*integral + renormalon)/r

        return total
    end

    function sigma_χ(; Q::Float64, χ::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String,
        μH_ratio::Float64=1.0, μJ_ratio::Float64=1.0, νJ_ratio::Float64=1.0, 
        μS_ratio::Float64=1.0, νS_ratio::Float64=1.0, bmax_ratio::Float64=1.0,
        parameters::Vector{Float64},
        Ω1::Float64, μ_ren_ratio::Float64=1.0)

        z = 0.5*(1-cos(χ))

        part = sigma_z(Q=Q, z=z, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL,
                        μH_ratio=μH_ratio, μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, 
                        μS_ratio=μS_ratio, νS_ratio=νS_ratio, bmax_ratio=bmax_ratio,
                        parameters=parameters,
                        Ω1=Ω1, μ_ren_ratio=μ_ren_ratio)

        total = 0.5*sin(χ)*part

        return total
    end

end