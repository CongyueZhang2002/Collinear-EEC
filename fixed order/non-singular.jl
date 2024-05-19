include("perturbation.jl")
include("SCET.jl")

function non_singular_func(; z::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    perturbation = perturbation_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)
    SCET = SCET_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    non_singular = perturbation - SCET

    return non_singular
end

function non_singular_sigma_χ(; χ::Float64, Q::Float64, μ_ini::Float64, αs_ini::Float64, order::Int64, nf::Int64)

    z = 0.5*(1-cos(χ))

    part = non_singular_func(z=z, Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, order=order, nf=nf)

    total = 0.5*sin(χ)*part

    return total
end

