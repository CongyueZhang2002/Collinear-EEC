using DataFrames

function ratios()

    ratio_matrix = DataFrame(μS=Float64[], μB=Float64[], νS=Float64[], νB=Float64[])

    for i in [1/2, 1, 2]
        for j in [1/2, 1, 2]
            for k in [1/2, 1, 2]
                for l in [1/2, 1, 2]    
                    push!(ratio_matrix, (μS=i, μB=j, νS=k, νB=l))
                end
            end
        end
    end

    n_rows = size(ratio_matrix, 1)

    constrained_matrix = DataFrame(μS=Float64[], μB=Float64[], νS=Float64[], νB=Float64[])

    for i in 1:n_rows

        μS = ratio_matrix[i, :μS]
        μB = ratio_matrix[i, :μB]
        νS = ratio_matrix[i, :νS]
        νB = ratio_matrix[i, :νB]

        if (
            (abs(μS/μB) < 3) && (abs(μS/νS) < 3) && (abs(νB/νS) < 3) &&
            (abs(μS/μB) > 1/3) && (abs(μS/νS) > 1/3) && (abs(νB/νS) > 1/3) &&
            (μS != 1 || μB != 1 || νS != 1 || νB != 1) && 
            (μB != 1/2 || νS != 1/2 || νB != 1) && (μB != 1/2 || νB != 2 || νS != 1) 
        )
            push!(constrained_matrix, (μS=μS, μB=μB, νS=νS, νB=νB))
        end
    end

    return constrained_matrix
end

display(ratios())