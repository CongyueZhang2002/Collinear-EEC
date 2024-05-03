using Distributed

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
    include("hard\\hard.jl")
    include("jet\\jet.jl")
    include("soft\\soft.jl") 
    include("sigma.jl")

end

function run_distributed_computation(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, μ::Float64, order::Int64, cores::Int64=4)
    try

        if length(workers()) < cores
            addprocs(cores - length(workers()))  
        end

        Q1 = Q
        μ_ini1 = μ_ini
        αs_ini1 = αs_ini
        nf1 = nf
        μ1 = μ
        order1 = order

        function compute_sigma_χ(χ)
            return sigma_χ(Q=Q1, χ=χ, μ_ini=μ_ini1, αs_ini=αs_ini1, nf=nf1, μ=μ1, order=order1)
        end

        results = pmap(compute_sigma_χ, χ_list)

        return results
    catch e
        println("An error occurred: ", e)
        rethrow(e)
    finally
        #rmprocs(workers()) 
    end
end
