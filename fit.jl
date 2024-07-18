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
include("sigma.jl")

function model(; xlist::Vector{Float64}, αs::Float64, Q::Float64, XLL::String, XLO::Int64,
                parameters::Vector{Float64}=[0.0,0.0,0.0], Ω1::Float64=0.0,
                μH_ratio::Float64=1.0, μJ_ratio::Float64=1.0, νJ_ratio::Float64=1.0, 
                μS_ratio::Float64=1.0, νS_ratio::Float64=1.0, bmax_ratio::Float64=1.0)

    χ_list = π / 180 * xlist

    resum = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=91.1876, αs_ini=αs, nf=5, XLL=XLL, 
                        parameters=parameters,Ω1=Ω1,                
                        μH_ratio=μH_ratio, μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, 
                        μS_ratio=μS_ratio, νS_ratio=νS_ratio, bmax_ratio=bmax_ratio)

    non_singular = Float64[]

    for x in χ_list
        push!(non_singular, non_singular_sigma_χ(χ=x, Q=Q, μ_ini=91.1876, αs_ini=αs, nf=5, order=XLO))
    end

    matched = resum + non_singular

    return matched
end