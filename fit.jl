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

function model(; xlist::Vector{Float64}, αs::Float64, parameters::Vector{Float64}=[0.0,0.0,0.0])

    χ_list = π / 180 * xlist

    resum = sigma_fast(Q=91.2, χ_list=χ_list, μ_ini=91.1876, αs_ini=αs, nf=5, order=2, parameters=parameters)

    non_singular = Float64[]

    for x in χ_list
        push!(non_singular, non_singular_sigma_χ(χ=x, Q=91.2, μ_ini=91.1876, αs_ini=αs, nf=5, order=2))
    end

    matched = resum + non_singular

    return matched
end