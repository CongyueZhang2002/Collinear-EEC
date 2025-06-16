function non_singular_func(; z::Float64, Q::Float64)
 
    perturbation = perturbation_func(z=z, Q=Q)

    FO = FO_func(z=z, Q=Q)

    non_singular = perturbation - FO
    
    return non_singular
end

function R_ratio(Q)
    αs = alpha_s_func(Q)
    r = 1 + αs/π + (αs/(4π))^2*(CF*CA*(123/2-44*z3) + CF*TF*nf*(-22+16*z3) - CF^2*3/2)
    return r
end

function non_singular_sigma_χ(; χ::Float64, Q::Float64)

    z = 0.5*(1-cos(χ))

    if z > 0.02
        global nloops_FO = nloops_NS_bulk 
    else
        global nloops_FO = nloops_NS_peak
    end

    part = non_singular_func(z=z, Q=Q)

    r = R_ratio(Q)

    total = 0.5*sin(χ)*part/r

    return total
end