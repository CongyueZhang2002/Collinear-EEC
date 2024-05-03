using Distributed
using HCubature
using SpecialFunctions

@everywhere include("square.jl")

function test(χ_list,b)

    if length(workers()) < 4
    addprocs(4 - length(workers()))  
    end

    b1=b

    function compute_sigma_χ(χ)
        return square(χ,b1)
    end

    results = pmap(compute_sigma_χ, χ_list)

    return results
end


