include("perturbation.jl")
include("SCET.jl")

function non_singular_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    if order < 3
        perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
        SCET = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

        non_singular = perturbation - SCET
    end

    if order == 3

        αs = alpha_s_func(μf=Q, μi=μ_ini, αs=αs_ini, order=4, nf=nf)

        if ((z<0.98) && (z>0.01))
            order3 =  7.5 * log(1-z)^5 
                    + 65 * log(1-z)^4 
                    + 204 * log(1-z)^3 
                    + 272 * log(1-z)^2 
                    + 154 * log(1-z) 
                    + 113 
                    + 0.3527 * log(z)^2/z 
                    - 7.747 * log(z)/z 
                    + 19.784/z
        else
            order3 = 0
        end

        perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=2, nf=nf)
        SCET = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=2, nf=nf)

        non_singular= perturbation - SCET + (αs/π)^3*order3
    end
    
    return non_singular
end

function non_singular_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = non_singular_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total
end

