function sigma_fast(; Q_list::Vector{Float64}, χ_list::Vector{Float64})

    @assert length(Q_list) == length(χ_list) "Q_list and χ_list must be the same length"
     
    function compute_sigma_χ(variables)
        Q, χ = variables
        return sigma_χ(Q=Q, χ=χ)
    end

    variables_list = zip(Q_list, χ_list)
    results = pmap(compute_sigma_χ, variables_list)

    return results

end

function model(; xlist::Vector{Float64}, Qlist::Vector{Float64})

    χ_list = π / 180 * xlist
    Q_list = Qlist

    resum = sigma_fast(Q_list=Q_list, χ_list=χ_list)

    non_singular = Float64[]

    for (Q, x) in zip(Q_list, χ_list)
        push!(non_singular, non_singular_sigma_χ(χ=x, Q=Q))
    end

    matched =  params[1]*resum + params[4]*non_singular

    return matched
end