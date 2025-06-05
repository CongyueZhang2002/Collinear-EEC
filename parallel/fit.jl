function sigma_fast(; Q_list::Vector{Float64}, χ_list::Vector{Float64}
    , μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
    params::Vector{Float64}=[0.0], μJ_ratio::Float64, μH_ratio::Float64)

    @assert length(Q_list) == length(χ_list) "Q_list and χ_list must be the same length"

    try
        μ_ini1 = μ_ini
        αs_ini1 = αs_ini
        nf1 = nf
        order1 = order
        μJ_ratio1 = μJ_ratio
        μH_ratio1 = μH_ratio        

        function compute_sigma_χ(variables)
            Q, χ = variables
            return sigma_χ(Q=Q, χ=χ, 
                            μ_ini=μ_ini1, αs_ini=αs_ini1, nf=nf1, order=order1,
                            params=params, μJ_ratio=μJ_ratio1, μH_ratio=μH_ratio1)
        end

        variables_list = zip(Q_list, χ_list)
        results = pmap(compute_sigma_χ, variables_list)

        return results
    catch e
        println("An error occurred: ", e)
        rethrow(e)
    finally
        #rmprocs(workers()) 
    end
end

function model(; xlist::Vector{Float64}, αs::Float64, Qlist::Vector{Float64}, order::Int64,
                params::Vector{Float64}=[0.0], μJ_ratio::Float64=1.0, μH_ratio::Float64=1.0)

    χ_list = π / 180 * xlist
    Q_list = Qlist

    resum = sigma_fast(Q_list=Q_list, χ_list=χ_list, μ_ini=91.1876, αs_ini=αs, nf=5, order=order, 
                        params=params, μJ_ratio=μJ_ratio, μH_ratio=μH_ratio )

    non_singular = Float64[]

    for (Q, x) in zip(Q_list, χ_list)
        push!(non_singular, non_singular_sigma_χ(χ=x, Q=Q, μ_ini=91.1876, αs_ini=αs, nf=5, order=order))
    end

    matched = resum + non_singular

    return matched
end