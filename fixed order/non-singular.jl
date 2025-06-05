function non_singular_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    #χ = acos(1-2*z)
    if order == 0    
        perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=1, nf=nf)
    elseif order == 1 
        perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=2, nf=nf)
    elseif order == 2
        perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=3, nf=nf)
    end

    SCET = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    non_singular = perturbation - SCET
    
    return non_singular
end

function R_ratio(; Q::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64)
    αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=2, nf=nf)
    r = 1 + αs/π + (αs/(4π))^2*(CF*CA*(123/2-44*z3) + CF*TF*nf*(-22+16*z3) - CF^2*3/2)
    return r
end

function non_singular_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    if z>0.02
        part = non_singular_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
    else
        if order < 2
            part = non_singular_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
        elseif order == 2
            part = non_singular_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=1, nf=nf)
        end
    end

    r = R_ratio(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf)

    total = 0.5*sin(χ)*part/r

    return total
end