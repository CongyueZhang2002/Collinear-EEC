using DataFrames
using Distributions
function custom_random()
    n = rand(Uniform(1.0,2.0)) 
    if rand() < 0.5
        return n  
    else
        return 1/n  
    end
end
function random_ratios(n)
    constrained_matrix = DataFrame(μS=Float64[], μB=Float64[], νS=Float64[], νB=Float64[])
    while size(constrained_matrix, 1) < n
        μS, μB, νS, νB = (custom_random(),custom_random(),custom_random(),custom_random())  # Generate random values

        # Check the constraints
        if (
            (abs(μS/μB) <= 2) && (abs(μS/νS) <= 2) && (abs(νB/νS) <= 2) &&
            (abs(μS/μB) >= 0.5) && (abs(μS/νS) >= 0.5) && (abs(νB/νS) >= 0.5)
        )
            push!(constrained_matrix, (μS=μS, μB=μB, νS=νS, νB=νB))
        end
    end

    return constrained_matrix
end
display(random_ratios(10))