# https://arxiv.org/pdf/1808.07867
using HCubature
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")

@everywhere function η_integral(; μi::Float64, μf::Float64, 
                      F0::Float64, F1::Float64, F2::Float64, F3::Float64, F4::Float64,
                      αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    β2 = β2_func(nf)   
    β3 = β3_func(nf)
    β4 = β4_func(nf)    

    αs_i = alpha_s_func(μf=μi, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
    αs_f = alpha_s_func(μf=μf, μi=μ_ini, αs=αs_ini, order=4,nf=nf)

    r = αs_f/αs_i
    B2 = β1^2/β0^2 - β2/β0
    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0
    B4 = β1^4/β0^4 - (3*β1^2*β2)/β0^3 + β2^2/β0^2 + (2*β1*β3)/β0^2 - β4/β0 

    η0 = - F0/(2β0)*log(r)

    η1 = - F0/(2β0)*αs_i/(4π)*(F1/F0 - β1/β0)*(r-1)

    η2 = - F0/(2β0)*(αs_i/(4π))^2*(B2 - F1*β1/(F0*β0) + F2/F0)*(r^2-1)/2  

    η3 = - F0/(2β0)*(αs_i/(4π))^3*(B3 + F1/F0*B2 - F2*β1/(F0*β0) + F3/F0)*(r^3-1)/3    

    η4 = - F0/(2β0)*(αs_i/(4π))^4*(B4 + F1/F0*B3 + F2/F0*B2 - F3*β1/(F0*β0) + F4/F0)*(r^4-1)/4

    if order == 1 
        total = η0
    elseif order == 2  
        total = η0 + η1
    elseif order == 3 
        total = η0 + η1 + η2
    elseif order == 4
        total = η0 + η1 + η2 + η3
    elseif order == 5
        total = η0 + η1 + η2 + η3 + η4
    end    

    return total
end

@everywhere function KΓ_integral(; μi::Float64, μf::Float64, scale::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    Q = scale # the paper's convention

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    β2 = β2_func(nf)   
    β3 = β3_func(nf)
    β4 = β4_func(nf)    

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)    

    αs_i = alpha_s_func(μf=μi, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
    αs_f = alpha_s_func(μf=μf, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
    αs_Q = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4,nf=nf)

    B2 = β1^2/β0^2 - β2/β0
    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0
    B4 = β1^4/β0^4 - (3*β1^2*β2)/β0^3 + β2^2/β0^2 + (2*β1*β3)/β0^2 - β4/β0

    function K0_func(α)
        K0 = Γ0/(4*β0^2)*(
            4*π*(α^(-1) + log(α)/αs_Q)
        )
        return K0
    end

    function K1_func(α)
        K1 = Γ0/(4*β0^2)*(
            (-2*π*αs_Q*β0*log(α)*(-2*β1*Γ0 + 2*β0*Γ1 + β1*Γ0*log(α/αs_Q^2)) + α*(-(β1*Γ0) + β0*Γ1)*(4*π*β0 + 2*αs_Q*β1 + αs_Q*β1*log(αs_Q/α)))/(4*π*αs_Q*β0^2*Γ0)
        )
        return K1
    end

    function K2_func(α)
        K2 = Γ0/(4*β0^2)*(
            (-4*π*αs_Q*(8*B2*π*αs_Q*β0^3*Γ0 + α^2*β1*(β1^2*Γ0 - β0*β2*Γ0 - β0*β1*Γ1 + β0^2*Γ2))*log(α) + α*(B2*αs_Q*β0*(32*π^2*β0^2*Γ0 + 4*π*(α - 2*αs_Q)*β0*(-(β1*Γ0) + β0*Γ1) + α*αs_Q*(-(β1^2*Γ0) + β0*β2*Γ0 + β0*β1*Γ1 - β0^2*Γ2)) + 2*π*(8*π*β0*(α*β1^2*Γ0 - α*β0*β2*Γ0 + 2*αs_Q*β0*β2*Γ0 - α*β0*β1*Γ1 + (α - 2*αs_Q)*β0^2*Γ2) + α*αs_Q*(3*β1^3*Γ0 - 3*β0*β1^2*Γ1 + 2*β0^2*β2*Γ1 + β0*β1*(-5*β2*Γ0 + 3*β0*Γ2))) + 4*π*α*αs_Q*β1*(β1^2*Γ0 - β0*β1*Γ1 + β0*(-(β2*Γ0) + β0*Γ2))*log(αs_Q)))/(128*π^3*αs_Q*β0^3*Γ0)
        )
        return K2
    end

    function K3_func(α)    
        K3 = Γ0/(4*β0^2)*(
            (α*(3*B3*αs_Q*β0*(96*π^3*α*β0^3*Γ0 + 16*π^2*(α^2 - 3*αs_Q^2)*β0^2*(-(β1*Γ0) + β0*Γ1) + 6*π*α*αs_Q^2*β0*(-(β1^2*Γ0) + β0*β1*Γ1 + β0*(β2*Γ0 - β0*Γ2)) + α^2*αs_Q^2*(β1^3*Γ0 - β0*β1^2*Γ1 + β0*β1*(-2*β2*Γ0 + β0*Γ2) + β0^2*(β3*Γ0 + β2*Γ1 - β0*Γ3))) + 8*π*α*(3*B2*α*αs_Q*β0*(4*π*β0*(β1^2*Γ0 - β0*β1*Γ1 + β0*(-(β2*Γ0) + β0*Γ2)) + αs_Q*(β1^3*Γ0 - β0*β1^2*Γ1 + β0*β1*(-2*β2*Γ0 + β0*Γ2) + β0^2*(β3*Γ0 + β2*Γ1 - β0*Γ3))) + 4*π*(α*αs_Q*(-4*β1^4*Γ0 + 4*β0*β1^3*Γ1 + β0*β1^2*(11*β2*Γ0 - 4*β0*Γ2) + 3*β0^2*(-(β2^2*Γ0) + β0*β3*Γ1 + β0*β2*Γ2) + β0^2*β1*(-7*β3*Γ0 - 7*β2*Γ1 + 4*β0*Γ3)) + 6*π*β0*(3*αs_Q*β0^2*(β3*Γ0 - β0*Γ3) - 2*α*(β1^3*Γ0 - β0*β1^2*Γ1 + β0*β1*(-2*β2*Γ0 + β0*Γ2) + β0^2*(β3*Γ0 + β2*Γ1 - β0*Γ3)))))) + 96*π^2*αs_Q*((-6*B3*π*αs_Q^2*β0^4*Γ0 + α^3*β1*(β1^3*Γ0 - β0*β1^2*Γ1 + β0*β1*(-2*β2*Γ0 + β0*Γ2) + β0^2*(β3*Γ0 + β2*Γ1 - β0*Γ3)))*log(α) - α^3*β1*(β1^3*Γ0 - β0*β1^2*Γ1 + β0*β1*(-2*β2*Γ0 + β0*Γ2) + β0^2*(β3*Γ0 + β2*Γ1 - β0*Γ3))*log(αs_Q)))/(18432*π^5*αs_Q*β0^4*Γ0)
        )
        return K3
    end

    function K4_func(α)
        K4 = Γ0/(4*β0^2)*(
            (B4*α*αs_Q*β0*(1024*π^4*α^2*β0^4*Γ0 + 192*π^3*(α^3 - 4*αs_Q^3)*β0^3*(-(β1*Γ0) + β0*Γ1) - 96*π^2*α*αs_Q^3*β0^2*(β1^2*Γ0 - β0*β1*Γ1 + β0*(-(β2*Γ0) + β0*Γ2)) - 16*π*α^2*αs_Q^3*β0*(-(β1^3*Γ0) + β0*β1^2*Γ1 + β0*β1*(2*β2*Γ0 - β0*Γ2) + β0^2*(-(β3*Γ0) - β2*Γ1 + β0*Γ3)) - 3*α^3*αs_Q^3*(β1^4*Γ0 - β0*β1^3*Γ1 + β0*β1^2*(-3*β2*Γ0 + β0*Γ2) + β0^2*β1*(2*β3*Γ0 + 2*β2*Γ1 - β0*Γ3) + β0^2*(β2^2*Γ0 - β0*(β4*Γ0 + β3*Γ1 + β2*Γ2) + β0^2*Γ4))) + 6*π*α^3*(3*B3*α*αs_Q*β0*(16*π^2*β0^2*(β1^2*Γ0 - β0*β1*Γ1 + β0*(-(β2*Γ0) + β0*Γ2)) + αs_Q^2*(-(β1^4*Γ0) + β0*β1^3*Γ1 + β0*β1^2*(3*β2*Γ0 - β0*Γ2) + β0^2*β1*(-2*β3*Γ0 - 2*β2*Γ1 + β0*Γ3) + β0^2*(-(β2^2*Γ0) + β0*β2*Γ2 + β0*(β4*Γ0 + β3*Γ1 - β0*Γ4)))) + 8*π*(3*B2*α*αs_Q*β0*(-4*π*β0*(β1^3*Γ0 - β0*β1^2*Γ1 + β0^2*(β3*Γ0 + β2*Γ1) + β0*β1*(-2*β2*Γ0 + β0*Γ2)) + 4*π*β0^4*Γ3 + αs_Q*(-(β1^4*Γ0) + β0*β1^3*Γ1 + β0*β1^2*(3*β2*Γ0 - β0*Γ2) + β0^2*β1*(-2*β3*Γ0 - 2*β2*Γ1 + β0*Γ3) + β0^2*(-(β2^2*Γ0) + β0*β2*Γ2 + β0*(β4*Γ0 + β3*Γ1 - β0*Γ4)))) + π*(3*α*αs_Q*(5*β1^5*Γ0 - 5*β0*β1^4*Γ1 + β0*β1^3*(-19*β2*Γ0 + 5*β0*Γ2) + β0^2*β1^2*(14*β3*Γ0 + 14*β2*Γ1 - 5*β0*Γ3) + 4*β0^3*(-2*β2*β3*Γ0 - β2^2*Γ1 + β0*β4*Γ1 + β0*β3*Γ2 + β0*β2*Γ3) + β0^2*β1*(13*β2^2*Γ0 - 9*β0*β2*Γ2 + β0*(-9*β4*Γ0 - 9*β3*Γ1 + 5*β0*Γ4))) + 16*π*β0*(4*αs_Q*β0^3*(β4*Γ0 - β0*Γ4) + 3*α*(β1^4*Γ0 - β0*β1^3*Γ1 + β0*β1^2*(-3*β2*Γ0 + β0*Γ2) + β0^2*β1*(2*β3*Γ0 + 2*β2*Γ1 - β0*Γ3) + β0^2*(β2^2*Γ0 - β0*(β4*Γ0 + β3*Γ1 + β2*Γ2) + β0^2*Γ4)))))) + 192*π^3*αs_Q*(-16*B4*π*αs_Q^3*β0^5*Γ0*log(α) - 3*α^4*β1*(β1^4*Γ0 - β0*β1^3*Γ1 + β0*β1^2*(-3*β2*Γ0 + β0*Γ2) + β0^2*β1*(2*β3*Γ0 + 2*β2*Γ1 - β0*Γ3) + β0^2*(β2^2*Γ0 - β0*(β4*Γ0 + β3*Γ1 + β2*Γ2) + β0^2*Γ4))*log(α/αs_Q)))/(589824*π^7*αs_Q*β0^5*Γ0)
        ) 
        return K4
    end

    L0 = K0_func(αs_f) - K0_func(αs_i)
    L1 = K1_func(αs_f) - K1_func(αs_i)
    L2 = K2_func(αs_f) - K2_func(αs_i)
    L3 = K3_func(αs_f) - K3_func(αs_i)
    L4 = K4_func(αs_f) - K4_func(αs_i)

    if order == 1
        total = L0
    elseif order == 2  
        total = L0 + L1
    elseif order == 3 
        total = L0 + L1 + L2
    elseif order == 4
        total = L0 + L1 + L2 + L3
    elseif order == 5
        total = L0 + L1 + L2 + L3 + L4
    end 

    return total
end

@everywhere function η_integral_pert(; μi::Float64, μf::Float64, 
                      F0::Float64, F1::Float64, F2::Float64, F3::Float64, F4::Float64,
                      αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    β2 = β2_func(nf)   
    β3 = β3_func(nf)
    β4 = β4_func(nf)    

    αs_i = alpha_s_func_pert(μf=μi, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
    αs_f = alpha_s_func_pert(μf=μf, μi=μ_ini, αs=αs_ini, order=4,nf=nf)

    r = αs_f/αs_i
    B2 = β1^2/β0^2 - β2/β0
    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0
    B4 = β1^4/β0^4 - (3*β1^2*β2)/β0^3 + β2^2/β0^2 + (2*β1*β3)/β0^2 - β4/β0 

    η0 = - F0/(2β0)*log(r)

    η1 = - F0/(2β0)*αs_i/(4π)*(F1/F0 - β1/β0)*(r-1)

    η2 = - F0/(2β0)*(αs_i/(4π))^2*(B2 - F1*β1/(F0*β0) + F2/F0)*(r^2-1)/2  

    η3 = - F0/(2β0)*(αs_i/(4π))^3*(B3 + F1/F0*B2 - F2*β1/(F0*β0) + F3/F0)*(r^3-1)/3    

    η4 = - F0/(2β0)*(αs_i/(4π))^4*(B4 + F1/F0*B3 + F2/F0*B2 - F3*β1/(F0*β0) + F4/F0)*(r^4-1)/4

    if order == 1 
        total = η0
    elseif order == 2  
        total = η0 + η1
    elseif order == 3 
        total = η0 + η1 + η2
    elseif order == 4
        total = η0 + η1 + η2 + η3
    elseif order == 5
        total = η0 + η1 + η2 + η3 + η4
    end    

    return total
end

#TEST

function KΓ_Integrand(; μ::Float64, Q::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs_μ = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    Γ = Γ_func(αs=αs_μ, order=order+1, nf=nf)
    
    return log(μ/Q)*Γ
end

function KΓ_test(; μi::Float64, μf::Float64, scale::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
    f(x)=1/x[1]*KΓ_Integrand(; μ=x[1], Q=scale, αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
    integral, error = hcubature(f, [μi], [μf])
    return integral
end

function η_Integrand(; μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs_μ = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

    Γ = Γ_func(αs=αs_μ, order=order+1, nf=nf)
    
    return Γ
end

function η_test(; μi::Float64, μf::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
    f(x)=1/x[1]*η_Integrand(μ=x[1], αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
    integral, error = hcubature(f, [μi], [μf])
    return integral
end

#@everywhere function KΓ_integral_old(; μi::Float64, μf::Float64, scale::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
#
#    Q = scale # the paper's convention
#
#    β0 = β0_func(nf)
#    β1 = β1_func(nf)   
#    β2 = β2_func(nf)   
#    β3 = β3_func(nf)
#    β4 = β4_func(nf)    
#
#    Γ0 = Γ0_func(nf)
#    Γ1 = Γ1_func(nf)
#    Γ2 = Γ2_func(nf)
#    Γ3 = Γ3_func(nf)
#    Γ4 = Γ4_func(nf)    
#
#    αs_i = alpha_s_func(μf=μi, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
#    αs_f = alpha_s_func(μf=μf, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
#
#    αs_Q = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4,nf=nf)
#
#    r = αs_f/αs_i
#    rQ = αs_i/αs_Q
#    B2 = β1^2/β0^2 - β2/β0
#    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0
#
#    K0 = Γ0/(4*β0^2) * (4π/αs_i) * (rQ*log(r) + 1/r - 1)
#
#    K1 = Γ0/(4*β0^2) * (
#          (Γ1/Γ0 - β1/β0) * (rQ*(r-1) - log(r))
#        - β1/(2β0) * (log(r)^2 + 2*log(rQ)*log(r))
#    )
#    
#    K2 = Γ0/(4*β0^2) * (αs_i/(4π)) * (
#          (Γ2/Γ0 - β1*Γ1/(β0*Γ0)) * ((1-r)^2/2 + (r^2-1)/2*(rQ-1))
#        + B2*(rQ*(r^2-1)/2 - log(r)/rQ) 
#        + (β1*Γ1/(β0*Γ0) - β1^2/β0^2) * ((r-1)*(1-log(rQ)) - r*log(r))
#    )
#    
#    K3 = Γ0/(4*β0^2) * ((αs_i/(4π))^2) * (
#        (Γ3/Γ0 - Γ2*β1/(Γ0*β0) + B2*Γ1/Γ0 + B3) * (r^3-1)/3*rQ
#        - β1/(2β0) * (Γ2/Γ0 - Γ1*β1/(Γ0*β0) + B2) * (r^2*log(r) + (r^2-1)*log(rQ))
#        - B3*log(r)/(2*rQ^2)
#        + (β3/β0 - β1*β2/β0^2 - 2*Γ3/Γ0 + 3*Γ2*β1/(Γ0*β0) - Γ1*β1^2/(Γ0*β0^2))*(r^2-1)/4 
#        + B2*(Γ1/Γ0 - β1/β0) * (1-r)/rQ
#    )
#    
#    if order == 0
#        total = K0
#    elseif order == 1  
#        total = K0 + K1
#    elseif order == 2 
#        total = K0 + K1 + K2
#    elseif order == 3
#        total = K0 + K1 + K2 + K3
#    end 
#
#    return total
#end