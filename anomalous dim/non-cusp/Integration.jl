# https://arxiv.org/pdf/1808.07867
using HCubature
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")

@everywhere function η_integral(; μi::Float64, μf::Float64, 
                      F0::Float64, F1::Float64, F2::Float64, F3::Float64,
                      αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    β2 = β2_func(nf)   
    β3 = β3_func(nf)    

    αs_i = alpha_s_func(μf=μi, μi=μ_ini, αs=αs_ini, order=order+1,nf=5)
    αs_f = alpha_s_func(μf=μf, μi=μ_ini, αs=αs_ini, order=order+1,nf=5)

    r = αs_f/αs_i
    B2 = β1^2/β0^2 - β2/β0
    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0

    η0 = - F0/(2β0)*log(r)

    η1 = - F0/(2β0)*αs_i/(4π)*(F1/F0 - β1/β0)*(r-1)

    η2 = - F0/(2β0)*(αs_i/(4π))^2*(B2 + F2/F0 - F1*β1/(F0*β0))*(r^2-1)/2  

    η3 = - F0/(2β0)*(αs_i/(4π))^3*(B3 + F1/F0*B2 - F2*β1/(F0*β0) + F3/F0)*(r^3-1)/3    

    if order == 0 
        total = η0
    end
    if order == 1  
        total = η0 + η1
    end
    if order == 2 
        total = η0 + η1 + η2
    end
    if order == 3
        total = η0 + η1 + η2 + η3
    end    

    return total
end

@everywhere function KΓ_integral(; μi::Float64, μf::Float64, scale::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    Q = scale # the paper's convention

    β0 = β0_func(nf)
    β1 = β1_func(nf)   
    β2 = β2_func(nf)   
    β3 = β3_func(nf)    

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)    

    αs_i = alpha_s_func(μf=μi, μi=μ_ini, αs=αs_ini, order=order+1,nf=5)
    αs_f = alpha_s_func(μf=μf, μi=μ_ini, αs=αs_ini, order=order+1,nf=5)

    αs_Q = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=order+1,nf=5)

    r = αs_f/αs_i
    rQ = αs_i/αs_Q
    B2 = β1^2/β0^2 - β2/β0
    B3 = -β1^3/β0^3 + 2*β1*β2/β0^2 - β3/β0

    K0 = Γ0/(4*β0^2) * (4π/αs_i) * (rQ*log(r) + 1/r - 1)

    K1 = Γ0/(4*β0^2) * (
          (Γ1/Γ0 - β1/β0) * (rQ*(r-1) - log(r))
        - β1/(2β0) * (log(r)^2 + 2*log(rQ)*log(r))
    )
    
    K2 = Γ0/(4*β0^2) * (αs_i/(4π)) * (
          (Γ2/Γ0 - β1*Γ1/(β0*Γ0)) * ((1-r)^2/2 + (r^2-1)/2*(rQ-1))
        + B2*(rQ*(r^2-1)/2 - log(r)/rQ) 
        + (β1*Γ1/(β0*Γ0) - β1^2/β0^2) * ((r-1)*(1-log(rQ)) - r*log(r))
    )
    
    K3 = Γ0/(4*β0^2) * ((αs_i/(4π))^2) * (
        (Γ3/Γ0 - Γ2*β1/(Γ0*β0) + B2*Γ1/Γ0 + B3) * (r^3-1)/3*rQ
        - β1/(2β0) * (Γ2/Γ0 - Γ1*β1/(Γ0*β0) + B2) * (r^2*log(r) + (r^2-1)*log(rQ))
        - B3*log(r)/(2*rQ^2)
        + (β3/β0 - β1*β2/β0^2 - 2*Γ3/Γ0 + 3*Γ2*β1/(Γ0*β0) - Γ1*β1^2/(Γ0*β0^2))*(r^2-1)/4 
        + B2*(Γ1/Γ0 - β1/β0) * (1-r)/rQ
    ) 

    if order == 0
        total = K0
    end
    if order == 1  
        total = K0 + K1
    end
    if order == 2 
        total = K0 + K1 + K2
    end
    if order == 3
        total = K0 + K1 + K2 + K3
    end 

    return total
end

#TEST

function KΓ_Integrand(; μ::Float64, Q::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs_μ = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

    Γ = Γ_func(αs=αs_μ, order=order+1, nf=nf)
    
    return log(μ/Q)*Γ
end

function KΓ_test(; μi::Float64, μf::Float64, scale::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
    f(x)=1/x[1]*KΓ_Integrand(; μ=x[1], Q=scale, αs_ini=0.118, μ_ini=91.2, order=3, nf=5)
    integral, error = hcubature(f, [μi], [μf])
    return integral
end

#print(" KΓ analytic: ",KΓ_integral(μi=10.0, μf=100.0, scale=0.01, αs_ini=0.118, μ_ini=91.2, order=3, nf=5))
#print(" KΓ numerical: ",KΓ_test(μi=10.0, μf=100.0, scale=0.01, αs_ini=0.118, μ_ini=91.2, order=3, nf=5))

function η_Integrand(; μ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs_μ = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

    Γ = Γ_func(αs=αs_μ, order=order+1, nf=nf)
    
    return Γ
end

function η_test(; μi::Float64, μf::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
    f(x)=1/x[1]*η_Integrand(μ=x[1], αs_ini=0.118, μ_ini=91.2, order=3, nf=5)
    integral, error = hcubature(f, [μi], [μf])
    return integral
end

a = η_integral(; μi=1., μf=100.0, 
F0=Γ0_func(5), F1=Γ1_func(5), F2=Γ2_func(5), F3=Γ3_func(5),
αs_ini=0.118, μ_ini=91.2, order=3, nf=5)

#print(" η analytic: ", a)
#print(" η numerical: ", η_test(μi=1.0, μf=100.0, αs_ini=0.118, μ_ini=91.2, order=3, nf=5))