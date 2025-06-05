function sigma_fast(; Q::Float64, χ_list::Vector{Float64}
    , μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
    params::Vector{Float64}=[0.0], Ω1::Float64=0.0)

    try
        Q1 = Q
        μ_ini1 = μ_ini
        αs_ini1 = αs_ini
        nf1 = nf
        order1 = order

        function compute_sigma_χ(χ)
            return sigma_χ(Q=Q1, χ=χ, 
                            μ_ini=μ_ini1, αs_ini=αs_ini1, nf=nf1, order=order1,
                            params=params)
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

function model(; xlist::Vector{Float64}, αs::Float64, Q::Float64, order::Int64,
                params::Vector{Float64}=[0.0])

    χ_list = π / 180 * xlist

    resum = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=91.1876, αs_ini=αs, nf=5, order=order, 
                        params=params)

    non_singular = Float64[]

    for x in χ_list
        push!(non_singular, non_singular_sigma_χ(χ=x, Q=Q, μ_ini=91.1876, αs_ini=αs, nf=5, order=order))
    end

    matched = resum + non_singular

    return matched
end